import networkx as nx
from revolutionhtl.nxTree import induced_colors
from src.polytomy_identification.TreePolytomies import TreePolytomies


def is_polytomi(T, node):
    return len(T[node]) > 2


def get_polytomies(T):
    return [node for node in T if is_polytomi(T, node)]


def get_trees_with_polytomies(gene_trees: nx.DiGraph) -> list[TreePolytomies]:
    trees_with_polytomies: list[TreePolytomies] = []

    for og, tree in enumerate(gene_trees):
        nodes_with_polytomies = get_polytomies(tree)

        if nodes_with_polytomies:
            X: dict[
                int, dict[
                    int, list[int | str]
                ]
            ] = {x: {} for x in nodes_with_polytomies}

            for x in nodes_with_polytomies:
                Y = list(tree.successors(x))
                for y_i in Y:
                    X[x][y_i] = induced_colors(tree, y_i, 'label')

            trees_with_polytomies.append(TreePolytomies(og, tree, X))

    return trees_with_polytomies
