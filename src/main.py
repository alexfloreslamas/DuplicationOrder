import csv
import pandas as pd
import networkx as nx
import Utils.Utils as utils
from src.Utils.Plots import plot
import src.neighbor_joining.DMSeries as dms
from revolutionhtl.nhxx_tools import get_nhx
import src.neighbor_joining.NanNeighborJoining as nnj


def extract_leaves_with_prefix(tree: nx.DiGraph) -> tuple[str, list[str]]:
    """
    Extracts the prefix from leaf names, checks if all leaves have the same prefix,
    and returns the suffixes of the leaves.

    Args:
        tree (nx.DiGraph): The input tree.

    Returns:
        tuple[str, list[str]]:
            - The common prefix (or None if prefixes mismatch).
            - The list of suffixes (empty list if prefixes mismatch).
    """
    leaves = [node for node in tree if tree.out_degree(node) == 0]
    leaf_labels = [tree.nodes[leaf].get('label', '') for leaf in leaves]
    prefixes = {label.split('|')[0] for label in leaf_labels}  # Extract prefixes

    if len(prefixes) == 1:  # All leaves have the same prefix
        prefix = prefixes.pop()
        suffixes = [label.split('|')[1] if '|' in label else '' for label in leaf_labels]  # Extract suffixes
        return prefix, suffixes
    return None, []


def computations(hits_path: str, trees_path: str, real_trees_base_path: str, output_file: str) -> None:
    distance_pairs, trees_with_polytomies = utils.load_distance_pairs_and_trees_with_polytomies(hits_path, trees_path)

    # TODO: Manually deleting the 58th tree since it's breaking the code. I'll check the causes tomorrow.
    del trees_with_polytomies[58]

    # Filtered structure to store trees and their leaves
    filtered_trees_with_polytomies = []

    # Filter trees by leaf prefix
    for tp in trees_with_polytomies:
        original_tree: nx.DiGraph = tp.get_tree()
        prefix, leaves = extract_leaves_with_prefix(original_tree)

        if prefix:  # Tree passes the filter
            filtered_trees_with_polytomies.append((tp, leaves))

    # Open the TSV file for writing
    with open(output_file, 'w', newline='') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t')
        writer.writerow([
            "og", "precision1", "recall1", "contradiction1", "precision2", "recall2", "contradiction2",
            "precision1==precision2", "recall1==recall2", "contradiction1==contradiction2"
        ])

        # Process filtered trees
        for tp, leaves in filtered_trees_with_polytomies:
            original_tree: nx.DiGraph = tp.get_tree()
            full_nx_resolved_tree: nx.DiGraph = original_tree.copy()

            # Find the corresponding real tree
            in_tree_newick: str = get_nhx(original_tree, name_attr='label')
            real_tree_file_name: str = f"g{utils.extract_file_name_from_newick(in_tree_newick)}.pruned.tree"
            re_tree_newick: str = utils.read_newick_from_file(real_trees_base_path, real_tree_file_name)
            real_tree = utils.custom_tree(re_tree_newick)

            # Compare leaves
            real_leaves = [node for node in real_tree if real_tree.out_degree(node) == 0]
            real_leaf_names = [real_tree.nodes[leaf].get('label', '') for leaf in real_leaves]

            if sorted(leaves) == sorted(real_leaf_names):  # Leaves match
                X: list[int] = tp.get_nodes_with_polytomies()

                for x in X:
                    Y: list[int] = tp.get_ys(x)
                    C: list[list[str]] = [tp.get_cluster(x, y_i) for y_i in Y]

                    # Compute the distance matrix for the NJ algorithm
                    D, _, info = dms.compute_distance_matrix(distance_pairs, C, Y)

                    if not utils.is_diagonal_zero_and_nan_elsewhere(D):
                        resolved_subtree_newick = nnj.resolve_tree_with_nan(
                            D, [str(y) if isinstance(y, int) else y for y in Y], x
                        )
                        full_nx_resolved_tree = utils.update_tree_with_newick(
                            full_nx_resolved_tree, node=x, newick_str=resolved_subtree_newick
                        )

                # Compute Newick and custom trees
                nj_tree_newick: str = utils.transform_newick(get_nhx(full_nx_resolved_tree, 1))
                in_custom_t = utils.custom_tree(in_tree_newick)
                nj_custom_t = utils.custom_tree(nj_tree_newick)
                re_custom_t = utils.custom_tree(re_tree_newick)

                # Compute metrics
                og = tp.get_og()
                precision1, recall1, contradiction1 = utils.get_precision_recall_contradiction(in_custom_t, re_custom_t)
                precision2, recall2, contradiction2 = utils.get_precision_recall_contradiction(nj_custom_t, re_custom_t)

                # Write results to the TSV file
                writer.writerow([
                    og, precision1, recall1, contradiction1, precision2, recall2, contradiction2,
                    precision1 == precision2, recall1 == recall2, contradiction1 == contradiction2
                ])

        print(f"For the output file: {output_file}, consider:")
        print("\t- precision1, recall1, contradiction1: Results of comparing (in_custom_t, re_custom_t)")
        print("\t- precision2, recall2, contradiction2: Results of comparing (nj_custom_t, re_custom_t)")


def main():
    # File paths
    hits_path:              str = '../input/tl_project_alignment_all_vs_all/'
    trees_path:             str = '../input/tl_project.reconciliation.tsv'
    real_trees_base_path:   str = "../input/true_gene_trees/"
    tsv_output_file:        str = "../output/results.tsv"                   # File to save the results
    plots_path:             str = "../output/plots/"                        # Path to save the plots

    #  -----------------------------------------------------------------------------------------------------------------

    computations(hits_path, trees_path, real_trees_base_path, tsv_output_file)
    df = pd.read_csv(tsv_output_file, sep='\t')
    plot(df, plots_path)


if __name__ == "__main__":
    main()
