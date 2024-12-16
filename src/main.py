import networkx as nx
import re
import Utils.Utils as utils
from pandas import read_csv
from revolutionhtl.nhxx_tools import read_nhxx, get_nhx
import src.neighbor_joining.DMSeries as dms
import src.neighbor_joining.NanNeighborJoining as nnj


if __name__ == "__main__":
    # Path to the hits file.
    hits_path = '../input/tl_project_alignment_all_vs_all/'             # Path to the hits file
    trees_path = '../input/tl_project.reconciliation.tsv'               # Path to trees

    distance_pairs = utils.load_hits_compute_distance_pairs(hits_path)  # Load distances
    gTrees = read_csv(trees_path, sep='\t')                             # Load trees
    gTrees = gTrees.set_index('OG').tree.apply(read_nhxx)               # Load trees

    trees_with_polytomies = utils.get_trees_with_polytomies(gTrees)     # Identify those trees with polytomies
    print(f"{len(trees_with_polytomies) = }")
    print(">>>"*11)

    # for idx, tp in enumerate(trees_with_polytomies):
    #     print(f"{idx = }")
    #     X: list[int] = tp.get_nodes_with_polytomies()
    #     print(tp)
    #     for x in X:
    #         Y: list[int] = tp.get_ys(x)
    #         C: list[list[str]] = [tp.get_cluster(x, y_i) for y_i in Y]
    #         D, info = dms.compute_distance_matrix(distance_pairs, C, Y)
    #         print(info)
    #
    #     print("---"*11)

    idx = 29
    print(f"{idx = }")
    tp = trees_with_polytomies[idx]
    print(tp)
    print("---"*11)

    X: list[int] = tp.get_nodes_with_polytomies()
    x: int = X[0]
    Y: list[int] = tp.get_ys(x)
    C: list[list[str]] = [tp.get_cluster(x, y_i) for y_i in Y]

    print(utils.info(X, x, Y, C))

    D, info = dms.compute_distance_matrix(distance_pairs, C, Y)
    print(info)

    print("---"*11)

    Y_str = [str(y) if isinstance(y, int) else y for y in Y]
    newick = nnj.resolve_tree_with_nan(D, Y_str, x)
    print(f"Newick format: {newick}\n")


    print(utils.tree_to_string(tp.get_tree(), 0, show_labels=True))

    solved_tree = utils.update_tree_with_newick(tp.get_tree(), node=x, newick_str=newick)
    print(utils.tree_to_string(solved_tree, 0, show_labels=True))


    print(get_nhx(solved_tree, 1))
