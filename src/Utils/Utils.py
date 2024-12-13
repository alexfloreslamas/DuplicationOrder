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
                    X[x][y_i] = induced_colors(tree, y_i, 'label')

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
