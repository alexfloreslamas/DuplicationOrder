import pandas
import numpy as np
from math import log
from revolutionhtl.parse_prt import load_all_hits_raw, normalize_scores


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


if __name__ == "__main__":
    # Path to the hits file.
    hits_path = '../input/tl_project_alignment_all_vs_all/'
    distance_pairs = load_hits_compute_distance_pairs(hits_path)

    # Example leaf identifiers
    l1: str = 'noD_5_3_1_1|G4_4'
    l2: str = 'noD_5_3_1_1|G1_1'
    l3: str = 'I do not exist'

    # Query distances
    print(get_pair_distance(distance_pairs, l1, l2))    # Expected: float, say abc.def...
    print(get_pair_distance(distance_pairs, l2, l1))    # Same as above, i.e, abc.def...
    print(get_pair_distance(distance_pairs, l1, l3))    # Expected: np.nan
    print(get_pair_distance(distance_pairs, l3, l1))    # Expected: np.nan
    print(get_pair_distance(distance_pairs, l2, l3))    # Expected: np.nan
    print(get_pair_distance(distance_pairs, l3, l2))    # Expected: np.nan
