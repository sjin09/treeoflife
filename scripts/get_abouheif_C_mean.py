#!/usr/bin/env python

"""
Calculate Abouheif's C_mean measure of phylogenetic signal

Notes:
    - Abouheif's matrix is a bistochastic matrix with non-null diagonal elements
        - All entries are non-negative
        - Each row sums to 1
        - Each column sums to 1
        - Here, Abouheif's matrix is a measure of phylogenetic proximity
        - Abouheif test measures the proximity between two taxa i and j
        as the inverse of the number of branches descending from each interior node
        in the path connecting i to j the root of the tree.
        - The proximity depends thus on the interior nodes in the path connecting i to j.
        - aii is a measure of how evolutionarily isolated species i is relative
        to other members of the phylogenetic tree under study.
        - The definition of A is independent of branch length and
        depends only on the topology of the tree.
    - Abouheif’s C mean tests for serial independence is based on the sum of the successive squared differences

Reference:
    - Pavoine, S., Ollier, S., Pontier, D., & Chessel, D. (2008). Testing for phylogenetic signal in phenotypic traits:
"""

import argparse
import csv
from collections import defaultdict
from pathlib import Path
import sys
from typing import Dict, List, Optional, Set, Tuple
import warnings

import ete3
import numpy as np
import pandas as pd
import plotnine as p9
from statsmodels.stats.multitest import multipletests

FILL_COLOURS = ["#F21D1D", "#2AA4BF"]
SAMPLE_LEVEL_TAXONOMIC_RANKS = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Sample"]
SPECIES_LEVEL_TAXONOMIC_RANKS = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]


def parse_args(args: Optional[List[str]] = None) -> argparse.Namespace:
    if args is None:
        args = sys.argv[1:]

    parser = argparse.ArgumentParser(
        description="Calculate Abouheif's C_mean measure of phylogenetic signal",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--taxonomic-classification",
        type=str,
        required=True,
        help="file to read taxonomic classification for each sample",
    )
    parser.add_argument(
        "--signature-exposures",
        type=str,
        required=True,
        help="file to read germline or somatic mutational signature exposures",
    )
    parser.add_argument("--rtol-sigs", type=str, required=True, help="file to read a list of artefactual signatures")
    parser.add_argument("--pdf", type=str, required=True, help="file to plot")
    parser.add_argument("-o", "--output", type=str, required=True, help="file to write")
    parser.add_argument("--is-somatic", required=False, action="store_true", help="is somatic mutational signature")
    args = args[1:]
    return parser.parse_args(args)


def get_branch_path(tree: ete3.TreeNode, node1: ete3.TreeNode, node2: ete3.TreeNode) -> List[ete3.TreeNode]:
    common_ancestor = tree.get_common_ancestor(node1, node2)
    path = [common_ancestor]
    path.extend(get_path_to_ancestor(common_ancestor, node1))
    path.extend(get_path_to_ancestor(common_ancestor, node2))

    while node1 != common_ancestor:
        node1 = node1.up
        path.append(node1)

    while node2 != common_ancestor:
        node2 = node2.up
        path.append(node2)

    return path


def get_exposure_per_signature_per_sample(
    signature_exposures_path: Path,
    rtol_sig_path: Path,
    is_somatic: bool,
) -> Tuple[Dict[str, Dict[str, float]], Set[str], List[str]]:
    # Initialize variable
    exposure_per_signature_per_sample = {}

    # Read the signature exposures file
    df = pd.read_csv(signature_exposures_path)

    # Load artefactual signatures
    rtol_sigs = [line.rstrip() for line in open(rtol_sig_path).readlines()]

    # Rename the columns
    if is_somatic:
        df.columns = ["Sample"] + ["sToL" + col for col in df.columns[1:]]
        df = df.drop(columns=rtol_sigs)  # Remove the rtol signatures
        col_names = ["Sample"] + [f"sToL{op_idx}" for op_idx, _col_name in enumerate(df.columns[1:])]
        with open("stol_signature_mapping_lookup.tsv", "w") as outfile:
            for (col_name_i, col_name_j) in zip(df.columns[1:], col_names[1:]):
                outfile.write(f"{col_name_i}\t{col_name_j}\n")
        df.columns = col_names
        df = df.drop(columns="sToL0")  # Remove the averaging component
    else:
        df.columns = ["Sample"] + ["gToL" + col for col in df.columns[1:]]
        df = df.drop(columns="gToL0")  # Remove the averaging component

    # Change sample names
    df["Sample"] = df["Sample"].apply(lambda x: x.split(".")[1] if "." in x else x)
    samples = set(df["Sample"])
    signatures = list(df.columns[1:])
    for (_idx, row) in df.iterrows():
        row_dict = row.drop("Sample").to_dict()
        exp_sum = sum(row_dict.values())
        normalised_row_dict = {}
        for (k, v) in row_dict.items():
            normalised_row_dict[k] = v / exp_sum
        exposure_per_signature_per_sample[row["Sample"]] = normalised_row_dict
    return exposure_per_signature_per_sample, samples, signatures


def get_path_to_ancestor(ancestor: ete3.TreeNode, node: ete3.TreeNode) -> List[ete3.TreeNode]:
    path = []
    while node != ancestor:
        node = node.up
        path.append(node)
    return path


def get_samples_per_species(taxonomic_classification_path: Path, samples: Set[str]) -> Dict[str, List[str]]:
    samples_per_species = defaultdict(list)
    df = pd.read_csv(taxonomic_classification_path)
    for (_idx, row) in df.iterrows():
        species = row["Species"]
        sample = row["Sample"]
        if sample not in samples:
            continue
        samples_per_species[species].append(sample)
    return dict(samples_per_species)


def get_sample_tree(taxonomic_classification_path: Path, samples: Set[str]) -> ete3.Tree:
    def get_child(parent_clade, child_name: str):
        """
        Find the child with the specified name or create a new one if not found.
        """
        for child in parent_clade.get_children():
            if child.name == child_name:
                return child
        new_child = parent_clade.add_child(name=child_name)
        return new_child

    # Initialize tree
    root = ete3.Tree(name="root")

    # Import and iterate through the DToL samplesheet
    df = pd.read_csv(taxonomic_classification_path)

    # Remove samples with incomplete taxonomic classification.
    # Check exactly the ranks used to build the tree -- not a positional column
    # slice, which is fragile to column reordering and previously also caught the
    # non-taxonomic "Common name" column, silently dropping ~25% of samples
    # (mostly plants, which have no common name) from the tree.
    df_subset = df[~(df[SAMPLE_LEVEL_TAXONOMIC_RANKS].eq(".").any(axis=1))]
    for (_idx, row) in df_subset.iterrows():
        sample = row["Sample"]
        # Don't add sample to the tree if the sample is absent from signature exposure spreadsheet
        if sample not in samples:
            continue
        current_node = root
        for taxonomic_rank in SAMPLE_LEVEL_TAXONOMIC_RANKS:
            sample_taxonomic_rank = row[taxonomic_rank]
            if sample_taxonomic_rank == ".":
                raise ValueError(f"Taxonomic rank {taxonomic_rank} is undefined for sample {sample}")
            current_node = get_child(current_node, sample_taxonomic_rank)
    return root


def get_species_tree(taxonomic_classification_path: Path, samples: Set[str]) -> ete3.Tree:
    def get_child(parent_clade, child_name: str):
        """
        Find the child with the specified name or create a new one if not found.
        """
        for child in parent_clade.get_children():
            if child.name == child_name:
                return child
        new_child = parent_clade.add_child(name=child_name)
        return new_child

    # Initialize tree
    root = ete3.Tree(name="root")

    # Import and iterate through the DToL samplesheet
    df = pd.read_csv(taxonomic_classification_path)

    # Remove samples with incomplete taxonomic classification.
    # As in get_sample_tree, check exactly the ranks used to build the tree.
    df_subset = df[~(df[SAMPLE_LEVEL_TAXONOMIC_RANKS].eq(".").any(axis=1))]
    for (_idx, row) in df_subset.iterrows():
        sample = row["Sample"]
        # Don't add sample to the tree if the sample is absent from signature exposure spreadsheet
        if sample not in samples:
            continue
        current_node = root
        for taxonomic_rank in SAMPLE_LEVEL_TAXONOMIC_RANKS:
            sample_taxonomic_rank = row[taxonomic_rank]
            current_node = get_child(current_node, sample_taxonomic_rank)
    return root


def get_tree(nwk_path: Path) -> ete3.Tree:
    return ete3.Tree(nwk_path)


def get_Abouheif_A_matrix(tree: ete3.TreeNode) -> np.ndarray:
    """
    Notes:
        - For a species i, aii is therefore the inverse of the product of
        the number of branches descending from each node from the species to the root.
        - For a couple (i, j), aij is the inverse of the product of
        the number of branches descending from each node in the path connecting i and j.
        - Abouheif proximity refers to the phylogenetic proximity underlying the test of Abouheif (see references).
             - Let P be the set of all the nodes in the path going from node1 to node2.
             - Let DDP be the number of direct descendants from each node in P.
             Then, the so-called ’Abouheif’ distance is the inverse of the product of all terms in DDP.
             - oriAbouheif returns a matrix with non-null diagonal elements, as formulated in Pavoine et al. (2008).
             - This matrix is bistochastic (all marginal sums equal 1),
             but this bilinear symmetric form does not give rise to a Moran’s index,
             since it requires a null diagonal.
            - Abouheif contains Abouheif’s proximities but has a null diagonal, giving rise to a Moran’s index.
    """
    # Get the leaves of the tree
    leaves = tree.get_leaves()
    n_leaves = len(leaves)

    # Initialize the matrix
    mtx = np.zeros((n_leaves, n_leaves))

    # Calculate the proximity scores for off-diagonal elements of the matrix
    for i in range(n_leaves):
        for j in range(n_leaves):
            if i != j:
                p = set(get_branch_path(tree, leaves[i], leaves[j]))
                ddp = [len(node.get_children()) for node in p]
                mtx[i, j] = 1 / np.prod(ddp)

    # Calculate the diagonal scores for originality of the species
    for k in range(n_leaves):
        mtx[k, k] = 1 - np.sum(mtx[k, :])  # row and columns sums to 1

    return mtx


def get_Abouheif_C_mean(
    tree: ete3.TreeNode,
    exp_per_sig_per_sample,
    sigs: List[str],
    samples_per_species: Dict[str, List[str]],
    iter: int = 1000,
) -> Tuple[List[float], List[float], List[float]]:
    """
    Notes:
        - Phylogenetic signal is the tendency of related species to resemble each other more
        than species drawn at random from the same tree.
         - Abouheif’s C mean tests for serial independence is based on the sum of the successive squared differences
         between trait values of neighbouring species (Abouheif 1999).
         - As there exist multiple ways to present the order of branches in a phylogenetic tree,
         Abouheif suggested Cmean as the mean value of a random subset of all possible representations.
         - Pavoine et al. (2008) provided an exact analytical value of the test.
         - They demonstrated that it uses Moran’s I statistic with a new matrix of phylogenetic proximities,
        which does not relate to branch length but focuses on topology and has a non-zero diagonal.

    C mean = (z^T A z / sum(A))

    Algorithm:
        - Calculate the signature exposure for each species
        - Perform Z-score normalisation
        - Calculate the score for phylogenetic signal
            - np.dot(traits, A_mtx) gives a vector of weighted trait values
            based on their phylogenetic distances.
            - np.dot(np.dot(traits, A_mtx, traits)) gives weighted sum of squared differences
            betwewen trait values, adjusted by phylogenetic distances.
        - Perform permutation test to calculate p values
        - Perform Benjamini-Hochberg correction for multiple hypothesis testingnp.do
    """

    # Set the seed for reproducibility
    np.random.seed(28)

    # Initialize variables
    C_means = []
    p_vals = []
    A_mtx = get_Abouheif_A_matrix(tree)

    # Calculate Abouheif's C_mean measure of phylogenetic signal for each signature
    for sig in sigs:
        # Get signature exposure for each species
        # species_sig_exps = []
        # for leaf_name in tree.iter_leaf_names():
        #     sample_sig_exps = []
        #     samples = samples_per_species[leaf_name]
        #     for sample in samples:
        #         sample_sig_exps.append(exp_per_sig_per_sample[sample][sig])
        #     species_sig_exps.append(np.mean(sample_sig_exps))

        # Perform Z-score normalisation
        # species_sig_exps = [
        #     (x - np.mean(species_sig_exps)) / np.std(species_sig_exps, ddof=0) for x in species_sig_exps
        # ]

        # Calculate the score for phylogenetic signal
        # C_mean = np.dot(np.dot(species_sig_exps, A_mtx), species_sig_exps) / np.sum(A_mtx)

        # Perform permutation test to calculate p values
        # null_C = np.zeros(iter)
        # for _ in range(iter):
        #     np.random.shuffle(species_sig_exps)
        #     null_C[_] = np.dot(np.dot(species_sig_exps, A_mtx), species_sig_exps) / np.sum(A_mtx)
        # p_val = np.sum(C_mean < null_C) / iter
        # C_means.append(C_mean)
        # p_vals.append(p_val)

        # Get signature exposure for each sample
        sample_sig_exps = [exp_per_sig_per_sample[leaf_name][sig] for leaf_name in tree.iter_leaf_names()]

        # Perform Z-score normalisation
        sample_sig_exps = [(x - np.mean(sample_sig_exps)) / np.std(sample_sig_exps, ddof=0) for x in sample_sig_exps]

        # Calculate the score for phylogenetic signal
        C_mean = np.dot(np.dot(sample_sig_exps, A_mtx), sample_sig_exps) / np.sum(A_mtx)

        # Perform permutation test to calculate p values
        null_C = np.zeros(iter)
        for _ in range(iter):
            np.random.shuffle(sample_sig_exps)
            null_C[_] = np.dot(np.dot(sample_sig_exps, A_mtx), sample_sig_exps) / np.sum(A_mtx)
        p_val = np.sum(C_mean < null_C) / iter
        C_means.append(C_mean)
        p_vals.append(p_val)

    # Perform Benjamini-Hochberg correction for multiple hypothesis testing
    _, q_vals, _, _ = multipletests(p_vals, method="fdr_bh")
    return C_means, p_vals, q_vals


def write_Abouheif_C_mean(
    sigs: List[str], C_means: List[float], p_vals: List[float], q_vals: List[float], output_path: Path
) -> None:
    fieldnames = ["Signature", "C_mean", "p_value", "q_value"]
    with open(output_path, "w") as outfile:
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        for (sig, C_mean, p_val, q_val) in zip(sigs, C_means, p_vals, q_vals):
            writer.writerow({"Signature": sig, "C_mean": C_mean, "p_value": p_val, "q_value": q_val})


def plot_Abouheif_C_mean(cmean_path: Path, pdf_path: Path, is_somatic: bool) -> None:
    df = pd.read_csv(cmean_path)
    df["q_bin"] = pd.cut(df["q_value"], bins=[0, 0.01, np.inf], labels=["< 0.01", ">= 0.01"], right=False)
    df["q_bin"] = df["q_bin"].cat.remove_unused_categories()

    # Reorder 'Signature' based on 'C_mean'
    df["Signature"] = pd.Categorical(
        df["Signature"],
        categories=df.groupby("Signature")["C_mean"].mean().sort_values(ascending=False).index,
        ordered=True,
    )

    # Generate base plot
    base_plot = (
        p9.ggplot(df, p9.aes(x="Signature", y="C_mean", fill="q_bin"))
        + p9.geom_bar(stat="identity")
        + p9.theme_classic(16)
        + p9.theme(axis_text_x=p9.element_text(rotation=90, hjust=1))
        + p9.guides(fill=None)
        + p9.scale_fill_manual(values=FILL_COLOURS)
    )

    if is_somatic:
        plot = base_plot + p9.labs(x="\nSomatic mutational signatures\n", y="\nAbouheif's C mean\n")
    else:
        plot = base_plot + p9.labs(x="\nGermline mutational signatures\n", y="\nAbouheif's C mean\n")

    # Generate the plot
    plot.save(pdf_path, width=22, height=12)


def get_phylogenetic_signal(
    taxonomic_classification_path: Path,
    signature_exposure_path: Path,
    rtol_sig_path: Path,
    pdf_path: Path,
    output_path: Path,
    is_somatic: bool,
):
    exp_per_sig_per_sample, samples, sigs = get_exposure_per_signature_per_sample(
        signature_exposure_path, rtol_sig_path, is_somatic
    )
    samples_per_species = get_samples_per_species(taxonomic_classification_path, samples)
    tree = get_sample_tree(taxonomic_classification_path, samples)
    C_means, p_vals, q_vals = get_Abouheif_C_mean(tree, exp_per_sig_per_sample, sigs, samples_per_species, 1000)
    write_Abouheif_C_mean(sigs, C_means, p_vals, q_vals, output_path)
    plot_Abouheif_C_mean(output_path, pdf_path, is_somatic)


def main() -> int:
    options = parse_args(sys.argv)
    get_phylogenetic_signal(
        taxonomic_classification_path=options.taxonomic_classification,
        signature_exposure_path=options.signature_exposures,
        rtol_sig_path=options.rtol_sigs,
        pdf_path=options.pdf,
        output_path=options.output,
        is_somatic=options.is_somatic,
    )
    return 0


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
