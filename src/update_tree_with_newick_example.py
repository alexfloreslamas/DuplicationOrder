import networkx as nx
import Utils.Utils as utils


if __name__ == "__main__":
    # Create the original DiGraph
    D = nx.DiGraph()
    D.add_edges_from([
        (0, 1), (1, 2), (2, 3), (2, 4), (1, 5),
        (5, 6), (6, 7), (6, 8), (8, 9), (8, 10),
        (5, 11), (11, 12), (11, 13), (13, 14), (13, 15),
        (5, 16), (16, 17), (16, 18), (18, 19), (18, 20)
    ])
    print(utils.tree_to_string(D, 0, show_labels=False))

    # Newick string representing the resolved subtree at node 5
    newick = "(16:7.845625481164477,(6:16.784747386860875,11:19.527174685311312):7.845625481164477)5;"


    # Resolve the polytomy and get a new graph
    new_D = utils.update_tree_with_newick(D, node=5, newick_str=newick)

    print(utils.tree_to_string(new_D, 0, show_labels=False))
