import numpy as np
import pandas as pd
import src.Utils.Utils as utils


def compute_distance_matrix(PD: pd.Series, C: list[list[str]], Y: list[str]) -> np.ndarray:
    """
    Compute the estimated distance matrix D based on the given "estimate" conditions.

    :param PD: pandas Series where the index is frozensets of IDs (tuples) and values are floats or np.nan.
    :param C: list of lists, where each sublist contains IDs corresponding to a cluster.
    :param Y: list of taxa labels corresponding to each cluster in C.
    :return: Symmetric distance matrix D as a 2D numpy array.
    """
    # Ensure Y matches the length of C
    if len(C) != len(Y):
        raise ValueError("The length of taxa labels (Y) must match the number of clusters (C).")

    k = len(C)  # Number of clusters
    D = np.zeros((k, k))  # Initialize the distance matrix with zeros

    # Iterate over all pairs of clusters (i, j)
    for i in range(k):
        for j in range(i, k):  # Compute only the upper triangle (i <= j)
            if i == j:
                D[i, j] = 0  # Diagonal is zero
            else:
                total = len(C[i]) * len(C[j])  # Initial total
                numerator = 0.0

                # Compute numerator and adjust total
                for z_i in C[i]:
                    for z_j in C[j]:
                        value = utils.get_pair_distance(PD, z_i, z_j)
                        if value is not None:
                            if not np.isnan(value):
                                numerator += value
                            else:  # PD[pair] = NaN
                                total -= 1
                        else:  # Pair not in PD
                            total -= 1

                # Avoid division by zero
                if total > 0:
                    D[i, j] = numerator / total
                else:
                    D[i, j] = np.nan  # Assign NaN if no valid pairs exist

                D[j, i] = D[i, j]  # Ensure symmetry

    return D


def test_compute_distance_matrix() -> None:
    # Example pairwise distances
    e0 = (
        # PD as pandas Series
        pd.Series({
            frozenset(("A", "B")): np.nan,
            frozenset(("A", "C")): np.nan,
            frozenset(("A", "D")): np.nan,
            frozenset(("B", "C")): np.nan,
            frozenset(("B", "D")): np.nan,
            frozenset(("C", "D")): np.nan
        }),
        # Clusters
        [
            ["A"],  # Cluster C1
            ["B"],  # Cluster C2
            ["C"],  # Cluster C3
            ["D"]   # Cluster C4
        ],
        # Example taxa labels
        ["y1", "y2", "y3", "y4"]
    )

    e1 = (
        # PD as pandas Series
        pd.Series({
            frozenset(("A", "B")): 2.0,
            frozenset(("A", "C")): 3.0,
            frozenset(("B", "C")): np.nan,
            frozenset(("A", "D")): 4.0,
            frozenset(("B", "D")): 5.0,
            frozenset(("X", "Y")): 10.0,  # Pair not related to the clusters
            frozenset(("Z", "A")): 1.5   # Pair not related to the clusters
        }),
        # Clusters
        [
            ["A", "B"],  # C1
            ["C"],       # C2
            ["D"]        # C3
        ],
        # Example taxa labels
        ["y1", "y2", "y3"]
    )

    e2 = (
        # PD as pandas Series
        pd.Series({
            frozenset(("A", "B")): 1.0,
            frozenset(("A", "C")): 2.5,
            frozenset(("A", "D")): 3.0,
            frozenset(("A", "E")): np.nan,
            frozenset(("B", "C")): 2.0,
            frozenset(("B", "D")): 2.8,
            frozenset(("B", "E")): 4.5,
            frozenset(("C", "D")): 1.5,
            frozenset(("C", "E")): np.nan,
            frozenset(("D", "E")): 3.5,
            frozenset(("F", "G")): 1.0,  # Additional pair unrelated to clusters
            frozenset(("H", "I")): np.nan  # Another unrelated pair
        }),
        # Clusters
        [
            ["A", "B", "C"],  # Cluster C1
            ["D", "E"],       # Cluster C2
            ["F"],            # Cluster C3 (single-element cluster)
        ],
        # Example taxa labels
        ["Cluster_1", "Cluster_2", "Cluster_3"]
    )

    e3 = (
        # PD as pandas Series
        pd.Series({
            frozenset(("A", "B")): 2.0,
            frozenset(("A", "C")): 3.0,
            frozenset(("A", "D")): 4.0,
            frozenset(("B", "C")): 5.0,
            frozenset(("B", "D")): 6.0,
            frozenset(("C", "D")): 7.0
        }),
        # Clusters
        [
            ["A", "B"],  # Cluster C1
            ["C"],       # Cluster C2
            ["D"]        # Cluster C3
        ],
        # Example taxa labels
        ["y1", "y2", "y3"]
    )

    test_cases = [
        e0,  # Fully disconnected clusters with all distances as np.nan
        e1,  # Mixed valid and missing distances with unrelated pairs in PD
        e2,  # Larger clusters with some missing distances
        e3   # Fully connected clusters with no missing distances
    ]

    # Compute distance matrix
    try:
        for PD, C, Y in test_cases:
            D = compute_distance_matrix(PD, C, Y)
            print("Taxa for matrix D:", Y)
            print("Distance matrix D:")
            print(f"{D}\n")
    except ValueError as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    test_compute_distance_matrix()