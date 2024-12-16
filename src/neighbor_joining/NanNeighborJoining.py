import numpy as np


def validate_input(D: np.ndarray) -> None:
    """
    Validate the input distance matrix.
    Ensures that the matrix is square and symmetric (allowing NaN).
    """
    if not (D.shape[0] == D.shape[1]):
        raise ValueError("Distance matrix must be square.")
    if not np.allclose(D, D.T, equal_nan=True):  # Allow symmetry with NaN
        raise ValueError("Distance matrix must be symmetric (allowing NaN).")


def identify_disconnected_nodes(all_D: np.ndarray, all_taxa: list[str]) -> tuple[list[str], list[str], np.ndarray]:
    """
    Identify disconnected nodes and filter the distance matrix for connected nodes.

    :param all_D: 2D numpy array of pairwise distances, possibly containing NaN.
    :param all_taxa: List of taxa names corresponding to the matrix rows/columns.
    :return: Tuple containing:
             - List of disconnected taxa,
             - List of connected taxa,
             - Filtered distance matrix for connected taxa.
    """
    disconnected_nodes = [
        all_taxa[i] for i in range(len(all_taxa)) if np.all(np.isnan(all_D[i, :][np.arange(len(all_taxa)) != i]))
    ]
    connected_mask = ~np.array([taxon in disconnected_nodes for taxon in all_taxa])
    connected_taxa = [taxon for taxon in all_taxa if taxon not in disconnected_nodes]
    filtered_D = all_D[connected_mask][:, connected_mask]
    return disconnected_nodes, connected_taxa, filtered_D


def neighbor_joining(D: np.ndarray, taxa: list[str]) -> dict[str, dict[str, float] | float]:
    """
    Perform the Neighbor Joining algorithm on the given distance matrix and taxa.

    :param D: 2D numpy array of pairwise distances (no NaN).
    :param taxa: List of taxa names corresponding to the matrix rows/columns.
    :return: Tree structure as a nested dictionary with branch lengths.
    """
    n: int = len(taxa)  # Number of taxa
    tree: dict[str, dict[str, float] | float] = {taxon: {} for taxon in taxa}

    while n > 2:
        # Step 1: Compute Q-matrix
        row_sums: np.ndarray = np.sum(D, axis=1)
        Q: np.ndarray = (n - 2) * D - row_sums[:, None] - row_sums

        # Step 2: Find the pair of distinct taxa i and j (i.e., i ≠ j) for which Q(i, j) is smallest
        np.fill_diagonal(Q, np.inf)
        i, j = np.unravel_index(np.argmin(Q), Q.shape)

        # Step 3: Calculate the distance from each of the taxa in the pair (i, j) to the new node u.
        u: str = f"({taxa[i]},{taxa[j]})"
        delta_i_u: float = 0.5 * D[i, j] + (row_sums[i] - row_sums[j]) / (2 * (n - 2))  # δ(i, u)
        delta_j_u: float = D[i, j] - delta_i_u                                          # δ(j, u)
        tree[u] = {taxa[i]: delta_i_u, taxa[j]: delta_j_u}

        # Step 4: Compute the distance from each taxon outside the pair (i, j) to the new node u,
        #         and update the D matrix to reflect this merge.
        new_distances: list[float] = [
            (D[i, k] + D[j, k] - D[i, j]) / 2 for k in range(len(D)) if k != i and k != j
        ]
        D = np.delete(D, (i, j), axis=0)
        D = np.delete(D, (i, j), axis=1)
        D = np.vstack((D, new_distances))
        new_distances.append(0)
        D = np.column_stack((D, new_distances))

        # Step 5: Update taxa list
        taxa.pop(max(i, j))
        taxa.pop(min(i, j))
        taxa.append(u)

        n -= 1

    # Final step: Add the last two clusters to the tree
    #             i.e., this is the root x
    tree[f"({taxa[0]},{taxa[1]})"] = {
        taxa[0]: D[0, 1] / 2,
        taxa[1]: D[0, 1] / 2,
    }

    return tree


def to_newick(tree_data: dict, root_name: str = "X") -> str:
    """
    Convert a nested dictionary tree to Newick format.

    :param tree_data: Dictionary with "tree" and "disconnected_nodes".
    :param root_name: Name of the root node to add.
    :return: Newick format string.
    """
    tree = tree_data["tree"]
    disconnected_nodes = tree_data["disconnected_nodes"]

    # Find the root key (the longest key)
    root = max(tree.keys(), key=len)

    # Recursive function to process nodes
    def recurse(node: str, subtree: dict, branch_length=None) -> str:
        if not subtree:  # Leaf node
            return f"{node}:{branch_length}" if branch_length is not None else node

        # Internal node: recursively process children
        children = []
        for child, length in subtree.items():
            children.append(recurse(child, tree.get(child, {}), length))
        combined = ",".join(children)

        # Add branch length if provided
        return f"({combined}){f':{branch_length}' if branch_length else ''}"

    # Start recursion from the root
    resolved_tree = recurse(root, tree[root])

    # Attach disconnected nodes to the root
    if disconnected_nodes:
        all_nodes = [resolved_tree] + [f"{node}:0" for node in disconnected_nodes]
        unrooted_newick = f"({','.join(all_nodes)})"
    else:
        unrooted_newick = resolved_tree

    # Add the specified root node
    # return f"({unrooted_newick}){root_name};"
    return f"{unrooted_newick}{root_name};"   # I believe this should be the new return


def resolve_tree_with_nan(full_D: np.ndarray, full_taxa: list[str], root_name: str) -> str:
    """
    Resolve a tree using Neighbor Joining, handling NaN values.

    :param full_D: 2D numpy array of pairwise distances, possibly containing NaN.
    :param full_taxa: List of taxa names corresponding to the matrix rows/columns.
    :param root_name: Name of the root node to add.
    :return: Newick format string.
    """
    # Validate input
    validate_input(full_D)

    # Handle the case where all distances are NaN
    if np.all(np.isnan(full_D)):
        # Create a Newick string with all taxa as leaves under the root_name
        leaves = ",".join([f"{taxon}:0" for taxon in full_taxa])
        return f"({leaves}){root_name};"

    # Identify disconnected nodes and filter the matrix
    disconnected_nodes, connected_taxa, filtered_D = identify_disconnected_nodes(full_D, full_taxa)

    filtered_D = filtered_D.astype(float)  # Convert to float if not already

    # Perform NJ on connected nodes
    tree = neighbor_joining(filtered_D, connected_taxa)
    print(f"Tree: {tree}")

    # Return tree and disconnected nodes
    return to_newick({"tree": tree, "disconnected_nodes": disconnected_nodes}, root_name)


def test_resolve_tree_with_nan() -> None:
    t0 = (
        ["A", "B", "C"],
        np.array([
            [np.nan, np.nan, np.nan],
            [np.nan, np.nan, np.nan],
            [np.nan, np.nan, np.nan],
        ])
    )

    t00 = (
        ["A", "B", "C", "D"],
        np.array([
            [np.nan, np.nan, np.nan, np.nan],
            [np.nan, np.nan, np.nan, np.nan],
            [np.nan, np.nan, np.nan, np.nan],
            [np.nan, np.nan, np.nan, np.nan],
        ])
    )

    t1 = (
        ["A", "B", "C"],
        np.array([
            [0, 5, np.nan],
            [5, 0, np.nan],
            [np.nan, np.nan, 0],
        ])
    )

    t2 = (
        ["A", "B", "C", "D"],
        np.array([
            [0, 5, np.nan, np.nan],
            [5, 0, np.nan, np.nan],
            [np.nan, np.nan, 0, np.nan],
            [np.nan, np.nan, np.nan, 0],
        ]),
    )

    t3 = (
        ["A", "B", "C", "D", "E"],
        np.array([
            [0, 5, 9, 9, 8],
            [5, 0, 10, 10, 9],
            [9, 10, 0, 8, 7],
            [9, 10, 8, 0, 3],
            [8, 9, 7, 3, 0],
        ])
    )

    test_cases: list[tuple[list[str], np.ndarray]] = [
        # Disconnected taxa
        t0,
        t00,
        # Partially connected taxa
        t1,
        t2,
        # Fully connected taxa
        t3
    ]

    root_node = "X"

    for input_taxa, input_D in test_cases:
        print(f"Newick format: {resolve_tree_with_nan(input_D, input_taxa, root_node)}\n")


if __name__ == "__main__":
    test_resolve_tree_with_nan()
