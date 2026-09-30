#!/usr/bin/env python

import argparse
from collections import defaultdict
import multiprocessing as mp
from pathlib import Path
import sys
from typing import Dict, List

import pysam


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="file to read expected SBS96 somatic mutation counts"
    )
    parser.add_argument(
        "--ref-fasta",
        type=Path,
        required=True,
        help="reference FASTA file to read"
    )
    parser.add_argument(
        "--target",
        type=Path,
        required=True,
        help="list of autosomes and sex chromosomes"
    )
    parser.add_argument(
        "--sample",
        type=str,
        required=True,
        help="sample id"
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=1,
        required=False,
        help="number of threads"
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=False,
        help="file to write mutation burden"
    )
    args = args[1:]
    return parser.parse_args(args)


NTS = ["A", "C", "G", "T"]
PUR = set(["A", "G"])
PUR_TO_PYR_LOOKUP = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
PYR = ["C", "T"]
TRINUCLEOTIDES = [f"{ntj}{nti}{ntk}" for nti in PYR for ntj in NTS for ntk in NTS]


def load_chroms(target_path: Path) -> List[str]:
    chroms = [line.strip() for line in open(target_path)]
    return chroms


def get_trinucleotide_count_per_chrom(
    chrom: str,
    sequence: str,
    ref_tri_count_per_tri_per_chrom: Dict[str, Dict[str, int]]
) -> None:
    ref_tri_count_per_tri = defaultdict(lambda: 0)
    for i in range(len(sequence)-2):
        base = sequence[i]
        if base == "N":
            continue
        trinucleotide = sequence[i:i+3]
        if trinucleotide[1] in PUR:
            trinucleotide_pyr = "".join(
                [PUR_TO_PYR_LOOKUP.get(base, "N") for base in trinucleotide[::-1]]
            )
            ref_tri_count_per_tri[trinucleotide_pyr] += 1
        else:
            ref_tri_count_per_tri[trinucleotide] += 1
    ref_tri_count_per_tri_per_chrom[chrom] = dict(ref_tri_count_per_tri)


def get_reference_sequence_trinucleotide_count(ref_fasta_path: Path, target_path, thread_count: int) -> Dict[str, int]:
    p = mp.Pool(thread_count)
    manager = mp.Manager()
    sequence_lookup = pysam.FastaFile(ref_fasta_path)
    if target_path:
        chroms = load_chroms(target_path)
    else:
        chroms = sequence_lookup.references
    ref_tri_count_per_tri_per_chrom = manager.dict()
    get_trinucleotide_count_per_chrom_arguments = [
        (
            chrom,
            sequence_lookup[chrom],
            ref_tri_count_per_tri_per_chrom
        )
        for chrom in chroms
    ]
    p.starmap(get_trinucleotide_count_per_chrom, get_trinucleotide_count_per_chrom_arguments)
    p.close()
    p.join()
    # get total counts
    ref_tri_count_per_tri = defaultdict(lambda: 0)
    for chrom in ref_tri_count_per_tri_per_chrom:
        for tri in TRINUCLEOTIDES:
            ref_tri_count_per_tri[tri] += ref_tri_count_per_tri_per_chrom[chrom].get(tri, 0)
    return dict(ref_tri_count_per_tri)


def load_sample_sbs96_counts(input_path: Path) -> Dict[str, int]:
    count_per_sbs96 = {}
    for line in open(input_path).readlines():
        if line.startswith("#"):
            continue
        elif line.startswith("sub"):
            continue
        fields = line.rstrip().split()
        sbs96 = fields[2]
        expected_count = float(fields[3])
        count_per_sbs96[sbs96] = expected_count
    return count_per_sbs96


def load_sample_trinucleotide_counts(input_path: Path) -> Dict[str, int]:
    count_per_tri = {}
    for line in open(input_path).readlines():
        if line.startswith("#"):
            continue
        elif line.startswith("sub"):
            continue
        fields = line.rstrip().split()
        tri = fields[1]
        sample_tri_count = int(fields[-1])
        count_per_tri[tri] = sample_tri_count
    return count_per_tri


def get_somatic_mutation_rate_per_trinucleotide(input_path: Path) -> Dict[str, float]:
    # get total number of somatic mutations per trinucleotide
    somatic_mutation_rate_per_tri = {}
    total_somatic_mutation_count_per_tri = defaultdict(lambda: 0)
    count_per_sbs96 = load_sample_sbs96_counts(input_path)
    sample_tri_count_per_tri = load_sample_trinucleotide_counts(input_path)
    for (sbs96, count) in count_per_sbs96.items():
        tri = f"{sbs96[0]}{sbs96[2]}{sbs96[6]}"
        total_somatic_mutation_count_per_tri[tri] += count
    # get somatic mutation rate per trinucleotide
    somatic_mutation_rate_per_tri = {}
    for tri in TRINUCLEOTIDES:
        mut_count = total_somatic_mutation_count_per_tri[tri]
        sample_tri_count = sample_tri_count_per_tri[tri]
        if sample_tri_count != 0:
            somatic_mutation_rate_per_tri[tri] = mut_count / sample_tri_count
        else:
            somatic_mutation_rate_per_tri[tri] = 0.0
    return somatic_mutation_rate_per_tri


def get_mutation_burden_per_cell(
    somatic_mutation_rate_per_tri: Dict[str, float],
    ref_tri_count_per_tri: Dict[str, int]
) -> float:
    burden = 0
    for tri in TRINUCLEOTIDES:
        burden += (somatic_mutation_rate_per_tri[tri] * ref_tri_count_per_tri[tri])
    return burden


def write_mutation_burden_per_cell(
    input_path: Path,
    ref_fasta_path: Path,
    target_path: Path,
    sample: str,
    thread_count: int,
    output_path: Path
) -> None:
    somatic_mutation_rate_per_tri = get_somatic_mutation_rate_per_trinucleotide(input_path)
    ref_tri_count_per_tri = get_reference_sequence_trinucleotide_count(
        ref_fasta_path,
        target_path,
        thread_count
    )
    burden_per_cell = get_mutation_burden_per_cell(somatic_mutation_rate_per_tri, ref_tri_count_per_tri)
    with open(output_path, "w") as outfile:
        outfile.write(f"{sample}\t{burden_per_cell}\n")


def main() -> int:
    options = parse_args(sys.argv)
    write_mutation_burden_per_cell(
        options.input,
        options.ref_fasta,
        options.target,
        options.sample,
        options.threads,
        options.output
    )
    return 0


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
