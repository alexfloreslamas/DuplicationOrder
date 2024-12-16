import pandas
import numpy as np
from math import log
import networkx as nx
from revolutionhtl.nxTree import induced_colors
from src.polytomy_identification.TreePolytomies import TreePolytomies
from revolutionhtl.parse_prt import load_all_hits_raw, normalize_scores


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
                    X[x][y_i] = list(induced_colors(tree, y_i, 'label'))

            trees_with_polytomies.append(TreePolytomies(og, tree, X))

    return trees_with_polytomies


def load_hits_compute_distance_pairs(hits_path: str) -> pandas.Series:
    """
    Load alignment hits, normalize scores, and compute pairwise distances.

    :param hits_path: Path to the hits file.
    :return: A pandas Series where the index is frozensets of leaf pairs and the values are distances.
    """
    df_hits = load_all_hits_raw(hits_path)  # Load alignment hits
    normalized_score = normalize_scores(df_hits, 'target')  # Compute distance
    distance = normalized_score.apply(
        lambda x: -log(min(x / 2, 1)) * 100
    )  # log correction of normalized bitscore a.k.a scoredist

    return distance


def get_pair_distance(distance_pairs_series: pandas.Series, leaf_1: str, leaf_2: str) -> float:
    """
    Retrieve the pairwise distance for two leaves from a pandas Series indexed by frozensets.

    :param distance_pairs_series: pandas Series where the index is frozensets of leaf pairs.
    :param leaf_1: The first leaf identifier.
    :param leaf_2: The second leaf identifier.
    :return: The distance between the pair of leaves if present, else np.nan.
    """
    pair = frozenset({leaf_1, leaf_2})
    return distance_pairs_series.get(pair, np.nan)


def info(X: list[int], x: int, Y: list[int], C: list[list[str]]) -> str:
    text: str = f"{X = }\n" \
                f"{x = }\n" \
                f"{Y = }\n" \
                f"C = [\n"
    for c in C:
        text += f"\t{c},\n"
    text += "]"
    return text


def tree_to_string(D: nx.DiGraph, node, level=0, prefix="", show_labels=True) -> str:
    # Get the label of the node if it exists, otherwise use the node ID
    label = D.nodes[node].get('label', str(node))

    # Create the current line of the tree
    connector = "|-" if level > 0 else ""  # Add '|-' only if not the root
    line = f"{prefix}{connector}{node} ({label})" if show_labels else f"{prefix}{connector}{node}"

    # Initialize the tree string with the current node
    tree_str = line + "\n"

    # Prepare the prefix for the next level
    new_prefix = prefix + ("| " if level > 0 else "  ")  # Add '| ' or '  '

    # Recurse for children and accumulate the result
    children = list(D.successors(node))
    for i, child in enumerate(children):
        is_last = i == len(children) - 1  # Adjust prefix for the last child
        child_prefix = new_prefix if not is_last else prefix + "  "
        tree_str += tree_to_string(D, child, level + 1, child_prefix, show_labels)

    return tree_str
