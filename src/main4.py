import csv
import networkx as nx
import Utils.Utils as utils
import src.neighbor_joining.DMSeries as dms
from revolutionhtl.nhxx_tools import get_nhx
import src.neighbor_joining.NanNeighborJoining as nnj


def main():
    # File paths
    hits_path = '../input/tl_project_alignment_all_vs_all/'
    trees_path = '../input/tl_project.reconciliation.tsv'
    real_trees_base_path = "../input/true_gene_trees/"
    output_file = "../output/results.tsv"  # File to save the results

    distance_pairs, trees_with_polytomies = utils.load_distance_pairs_and_trees_with_polytomies(hits_path, trees_path)

    # TODO: Manually deleting the 58th tree since it's breaking the code. I'll check the causes tomorrow.
    del trees_with_polytomies[58]

    # Open the TSV file for writing
    with open(output_file, 'w', newline='') as tsvfile:
        # Create a TSV writer
        writer = csv.writer(tsvfile, delimiter='\t')

        # Write the header row
        writer.writerow(["og", "precision1", "recall1", "contradiction1", "precision2", "recall2", "contradiction2"])

        # Process each tree with polytomies
        for idx, tp in enumerate(trees_with_polytomies):
            # print(f"{idx = }")
            original_tree: nx.DiGraph = tp.get_tree()
            full_nx_resolved_tree: nx.DiGraph = original_tree.copy()

            X: list[int] = tp.get_nodes_with_polytomies()

            for x in X:
                Y: list[int] = tp.get_ys(x)
                C: list[list[str]] = [tp.get_cluster(x, y_i) for y_i in Y]

                # Computes the distance matrix for the Neighbor-Joining (NJ) algorithm.
                D, _, info = dms.compute_distance_matrix(distance_pairs, C, Y)

                if not utils.is_diagonal_zero_and_nan_elsewhere(D):
                    resolved_subtree_newick = nnj.resolve_tree_with_nan(D, [str(y) if isinstance(y, int) else y for y in Y], x)
                    full_nx_resolved_tree = utils.update_tree_with_newick(
                        full_nx_resolved_tree, node=x, newick_str=resolved_subtree_newick
                    )

            # Compute Newick strings and custom trees
            in_tree_newick: str = get_nhx(original_tree, name_attr='label')
            nj_tree_newick: str = utils.transform_newick(get_nhx(full_nx_resolved_tree, 1))
            real_tree_file_name: str = f"g{utils.extract_file_name_from_newick(in_tree_newick)}.pruned.tree"
            re_tree_newick: str = utils.read_newick_from_file(real_trees_base_path, real_tree_file_name)

            in_custom_t = utils.custom_tree(in_tree_newick)
            nj_custom_t = utils.custom_tree(nj_tree_newick)
            re_custom_t = utils.custom_tree(re_tree_newick)

            # Compute metrics
            og = tp.get_og()
            precision1, recall1, contradiction1 = utils.get_precision_recall_contradiction(in_custom_t, re_custom_t)
            precision2, recall2, contradiction2 = utils.get_precision_recall_contradiction(nj_custom_t, re_custom_t)

            # print(f"\t{og = }")
            # print(f"\tPerformance(in_t, re_t): {precision1, recall1, contradiction1}")
            # print(f"\tPerformance(nj_t, re_t): {precision2, recall2, contradiction2}")

            # Write the metrics to the TSV file
            writer.writerow([og, precision1, recall1, contradiction1, precision2, recall2, contradiction2])

        print(f"For the output file: {output_file}, consider:")
        print("\t- precision1, recall1, contradiction1: Results of comparing (in_custom_t, re_custom_t)")
        print("\t- precision2, recall2, contradiction2: Results of comparing (nj_custom_t, re_custom_t)")


if __name__ == "__main__":
    main()
