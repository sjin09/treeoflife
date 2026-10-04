#!/usr/bin/env python

import argparse
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional
import sys

from ete3 import Tree
import pandas as pd

TAXONOMIC_RANKS = ["Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]


def parse_args(args: Optional[List[str]] = None) -> argparse.Namespace:
    if args is None:
        args = sys.argv[1:]

    parser = argparse.ArgumentParser(
        description="Get newick tree from taxonomic classification",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="file to read taxonomic classifications per sample",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="file to write newick tree (.nwk)",
    )
    return parser.parse_args(args)


def get_sample_count_per_reference_sample(df: pd.DataFrame) -> Dict[str, int]:
    # Count number of samples for each reference sample
    sample_count_per_reference_sample = defaultdict(lambda: 0)
    for (_idx, row) in df.iterrows():
        ref_sample = row["Reference sample"]
        sample_count_per_reference_sample[ref_sample] += 1
    return sample_count_per_reference_sample 


def get_sample_count_per_species(df: pd.DataFrame) -> Dict[str, int]:
    # Count number of samples for each reference sample
    sample_count_per_species = defaultdict(lambda: 0)
    for (_idx, row) in df.iterrows():
        species_name = row["Species"]
        sample_count_per_species[species_name] += 1
    return sample_count_per_species

# def quote_newick_label(s: str) -> str:
#     # Newick escapes single quotes by doubling them
#     # s = s.replace("'", "''")
#     return f"'{s}'"

def write_newick_tree(input_path: Path, output_path: Path) -> None:
    def get_child(parent_clade, child_name: str):
        """
        Find the child with the specified name or create a new one if not found.
        """
        for child in parent_clade.get_children():
            if child.name == child_name:
                return child
        new_child = parent_clade.add_child(name=child_name)
        return new_child
    
    # Import and iterate through the DToL samplesheet
    df = pd.read_csv(input_path)

    # Subset data frame # Remove samples with incomplete taxonomic classification
    df_subset = df[~(df.iloc[:, 3:10].eq(".").any(axis=1))]

    # Calculate number of sample for each reference sample
    sample_count_per_species = get_sample_count_per_species(df_subset)
    sample_count_per_reference_sample = get_sample_count_per_reference_sample(df_subset)

    # Initialize tree
    root = Tree(name="root")

    # Build newick tree
    for (_idx, row) in df_subset.iterrows():
        current_node = root
        sample = row["Sample"]
        species = row["Species"]
        ref_sample = row["Reference sample"]
        sample_taxonomic_ranks = list(row[TAXONOMIC_RANKS])
        for taxonomic_rank, sample_taxonomic_rank in zip(TAXONOMIC_RANKS, sample_taxonomic_ranks):
            if taxonomic_rank == "Species":
                if ref_sample == sample:
                    if sample_count_per_reference_sample[ref_sample] == 1:
                        if sample_count_per_species[species] == 1:
                            current_node = get_child(current_node, sample_taxonomic_rank)
                        else:
                            sample_taxonomic_rank = f"{sample_taxonomic_rank} ({sample})"
                            current_node = get_child(current_node, sample_taxonomic_rank)
                    else:
                        sample_taxonomic_rank = f"{sample_taxonomic_rank} ({sample})"
                        current_node = get_child(current_node, sample_taxonomic_rank)
                else:
                    sample_taxonomic_rank = f"{sample_taxonomic_rank} ({sample})"
                    current_node = get_child(current_node, sample_taxonomic_rank)
            else:
                current_node = get_child(current_node, sample_taxonomic_rank)

    # Export tree to newick format
    root.write(format=2, outfile=str(output_path))
    return root


def plot_circular_tree(tree: Tree, output_path: Path, leaf_fontsize: int = 7, canvas_px: int = 2800) :
    from ete3 import TreeStyle, NodeStyle, TextFace  # needs PyQt5; imported here so the newick is written without it

    # Apply default thin branches everywhere first
    base_ns = NodeStyle()  # Create a default node style object
    base_ns["fgcolor"] = "#333333"  # Set the default branch/line color to dark gray
    base_ns["hz_line_width"] = 0.8  # Set horizontal branch line width
    base_ns["vt_line_width"] = 0.8  # Set vertical branch line width (thickness)
    base_ns["size"] = 0   # Remove node symbols (no dots at nodes)
    for n in tree.traverse():  # Traverse all nodes in the tree
        n.set_style(base_ns)  # Apply the default style to each node

    # TreeStyle for circular layout
    ts = TreeStyle()  # Create a TreeStyle object to control the rendering
    ts.mode = "c"  # Use circular ("c") layout mode
    # ts.show_leaf_name = Tree  # Do not automatically show leaf names 
    ts.show_leaf_name = False  # Do not automatically show leaf names 
    ts.arc_start = -180  # Start the circle drawing arc at -180 degrees
    ts.arc_span = 360  # Span the arc over a full circle (360 degrees)
    ts.optimal_scale_level = "full"   # Scale the tree so it fully fits into the drawing canvas
    ts.branch_vertical_margin = 8  # Set vertical spacing between branches (affects readability in circular trees)

    # Add species labels to leaves
    for leaf in tree.iter_leaves():
        # species name is leaf.name since we built path down to Species
        face = TextFace(leaf.name, fsize=leaf_fontsize, fstyle="italic")
        # position the label at the end of the branch (aligned radially)
        leaf.add_face(face, column=0, position="aligned")
    
    tree.render(output_path, w=canvas_px, units="px", tree_style=ts)

def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input CSV")
    parser.add_argument("-o", "--output", required=True, help="Output Newick")
    options = parse_args()
    tree = write_newick_tree(options.input, options.output)
    try:
        plot_circular_tree(tree, "{}.pdf".format(options.output.stem))
    except ImportError as e:
        print(f"Skipped circular tree PDF: {e}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
