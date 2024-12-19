import Utils.Utils as utils
from pandas import read_csv
import src.neighbor_joining.DMSeries as dms
from revolutionhtl.nhxx_tools import read_nhxx, get_nhx
import src.neighbor_joining.NanNeighborJoining as nnj


def main():
    # File paths
    hits_path = '../../input/tl_project_alignment_all_vs_all/'
    trees_path = '../../input/tl_project.reconciliation.tsv'
    real_trees_base_path = "../../input/true_gene_trees/"

    # Load the hits and gtrees data from input files
    distance_pairs = utils.load_hits_compute_distance_pairs(hits_path)  # Load distances
    gTrees = read_csv(trees_path, sep='\t')                             # Load trees
    gTrees = gTrees.set_index('OG').tree.apply(read_nhxx)               # Load trees

    print(f"{len(gTrees) = }")
    trees_with_polytomies = utils.get_trees_with_polytomies(gTrees)     # Identify those trees with polytomies
    print(f"{len(trees_with_polytomies) = }\n\n{'-'*80}\n")

    # utils.print_all_trees_with_polytomies(trees_with_polytomies, distance_pairs, verbose=True)

    # Analyze the polytomy at a specific index
    idx = 29
    # idx = 32
    # idx = 55
    print(f"{idx = }")
    tp = trees_with_polytomies[idx]
    print(f"{tp}\n{'-'*80}\n")

    X: list[int] = tp.get_nodes_with_polytomies()
    x: int = X[0]
    Y: list[int] = tp.get_ys(x)
    C: list[list[str]] = [tp.get_cluster(x, y_i) for y_i in Y]
    print(utils.info(X, x, Y, C))

    # Computes the distance matrix for the Neighbor-Joining (NJ) algorithm.
    # TODO: Handle the case where the polytomy forms a simple star tree.
    # In such cases, distances should be taken as-is (raw distances), and no "distance estimates" should be computed.
    # This can likely be addressed with a conditional check (e.g., an 'if' statement).
    D, _, info = dms.compute_distance_matrix(distance_pairs, C, Y)
    print(f"{info}\n{'-'*80}\n")

    # Resolve polytomy using NJ
    Y_str = [str(y) if isinstance(y, int) else y for y in Y]
    resolved_subtree_newick = nnj.resolve_tree_with_nan(D, Y_str, x)
    print(f"{resolved_subtree_newick = }\n")

    # Print in-tree and resolved out-tree
    rho: int = 0
    original_tree = tp.get_tree()
    print(f"in-tree:\n{utils.tree_to_string(original_tree, rho, show_labels=True)}")        # full tree with polytomies

    # Creates a new nx.DiGraph tree with the polytomy at node 'x' resolved using the provided Newick string.
    # The original tree remains unchanged.
    full_nx_resolved_tree = utils.update_tree_with_newick(original_tree, node=x, newick_str=resolved_subtree_newick)
    print(f"nj-tree:\n{utils.tree_to_string(full_nx_resolved_tree, rho, show_labels=True)}")    # full resolved tree

    # Newick representations
    in_tree_newick:         str = get_nhx(original_tree, name_attr='label')
    nj_tree_newick:         str = utils.transform_newick(get_nhx(full_nx_resolved_tree, 1))
    real_tree_file_name:    str = f"g{utils.extract_file_name_from_newick(in_tree_newick)}.pruned.tree"
    re_tree_newick:         str = utils.read_newick_from_file(real_trees_base_path, real_tree_file_name)

    print(f"in-Newick: {in_tree_newick}")                                                       # tree with polytomies
    print(f"nj-Newick: {nj_tree_newick}")                                                       # resolved tree
    print(f"re-Newick: {re_tree_newick}")                                                       # real tree


if __name__ == "__main__":
    main()
