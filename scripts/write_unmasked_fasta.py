#!/usr/bin/env python

import argparse
from pathlib import Path
import sys

import pysam


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=doc, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="FASTA file to read"
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="file to write"
    )
    args = args[1:]
    return parser.parse_args(args)


def write_unmasked_fasta(input_path: Path, output_path: Path):
    sequence_lookup = pysam.FastaFile(input_path)
    with open(output_path, "w") as out_fh:
        for chrom in sequence_lookup.references:
            sequence = sequence_lookup.fetch(chrom)
            unmasked = sequence.upper()
            out_fh.write(f">{chrom}\n")

            # Write sequence in 60-char lines (standard FASTA)
            for i in range(0, len(unmasked), 60):
                out_fh.write(unmasked[i:i+60] + "\n")
    sequence_lookup.close()

def main() -> int:
    options = parse_args(sys.argv)
    write_unmasked_fasta(options.input, options.output)
    return 0


if name == "main":
    exit_code = main()
    sys.exit(exit_code)