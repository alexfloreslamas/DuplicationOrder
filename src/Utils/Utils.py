import os
import re
import math
import pandas
import numpy as np
from math import log
import networkx as nx
from Bio import Phylo
from io import StringIO
from pandas import read_csv
import src.neighbor_joining.DMSeries as dms
from revolutionhtl.nhxx_tools import read_nhxx
from revolutionhtl.nxTree import induced_colors
from itertools import chain, product, combinations
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


def update_tree_with_newick(D, node, newick_str) -> nx.DiGraph:
    """
    Resolve a polytomy in a NetworkX DiGraph using a Newick string and return a new graph.

    Args:
        D (nx.DiGraph): Original directed graph.
        node (int): Node with a polytomy to resolve.
        newick_str (str): Newick string representing the resolved subtree.

    Returns:
        nx.DiGraph: A new graph with the resolved polytomy.
    """
    # Create a copy of the original graph
    new_graph = D.copy()

    # Parse the Newick string into a tree
    handle = StringIO(newick_str)
    tree = Phylo.read(handle, "newick")

    # Generate new nodes for internal nodes introduced in the Newick tree
    new_node_id = max(new_graph.nodes) + 1  # Start creating new nodes from max existing ID + 1
    new_node_map = {}

    def map_newick_clade(clade):
        nonlocal new_node_id
        if clade.is_terminal():
            # Terminal nodes are the same as the original nodes
            return int(clade.name)
        else:
            # Create a new internal node
            new_internal_node = new_node_id
            new_node_id += 1
            new_node_map[clade] = new_internal_node
            return new_internal_node

    # Traverse the Newick tree to add the resolved structure
    def add_newick_edges(clade, parent):
        if clade.is_terminal():
            # Terminal node: add edge to the terminal node
            return int(clade.name)
        else:
            # Internal node: add edge to a new internal node
            internal_node = map_newick_clade(clade)
            for child in clade.clades:
                child_node = add_newick_edges(child, internal_node)
                new_graph.add_edge(internal_node, child_node)
            return internal_node

    # Remove old edges from the polytomy node in the new graph
    for child in list(new_graph.successors(node)):
        new_graph.remove_edge(node, child)

    # Add the resolved structure to the graph
    root_clade = tree.clade
    clades = root_clade.clades

    if len(clades) >= 1:
        # First child of the Newick tree becomes a direct child of `node`
        first_child = add_newick_edges(clades[0], node)
        new_graph.add_edge(node, first_child)

    if len(clades) > 1:
        # Second clade introduces a new internal node
        new_internal_node = new_node_id
        new_node_id += 1
        new_graph.add_edge(node, new_internal_node)
        # Add the 'label' attribute to the new_internal_node.
        new_graph.nodes[new_internal_node]['label'] = 'D'       # Discussed with Toño. 17/12/2024
        new_graph.nodes[new_internal_node]['node_id'] = None    # TODO: Could be new_internal_node - 1
        new_graph.nodes[new_internal_node]['species'] = None

        # Add the subtrees rooted at the second clade
        for child_clade in clades[1].clades:
            child_node = add_newick_edges(child_clade, new_internal_node)
            new_graph.add_edge(new_internal_node, child_node)

    return new_graph

def transform_newick(input_newick):
    """
    Transforms a Newick string by reordering attributes inside square brackets.

    Specifically, the function:
      - Reorders attributes in the format: 'label[node_id=...;species=...]'
      - Handles attributes appearing in any order: species, node_id, label
      - Ensures the output format always starts with 'label' followed by 'node_id' and 'species'.

    Args:
        input_newick (str): The original Newick string.

    Returns:
        str: The transformed Newick string with reordered attributes.
    """

    # Function to parse and reorder attributes
    def reorder_attributes(match):
        """
        Extracts and reorders attributes found within square brackets.

        Args:
            match (re.Match): A regex match object containing attributes inside square brackets.

        Returns:
            str: A string with reordered attributes in the desired format.
        """
        content = match.group(1)  # Get the full content inside the brackets
        # Split the attributes by ';' and create a key-value dictionary
        attributes = dict(kv.split('=') for kv in content.split(';') if kv)

        # Extract specific attributes (order is important)
        species = attributes.get('species', None)
        node_id = attributes.get('node_id', None)
        label = attributes.get('label', None)

        # Construct the replacement string in the desired order
        replacement = f"{label or ''}[node_id={node_id or ''}"
        if species:
            replacement += f";species={species}"
        replacement += "]"

        return replacement

    # Regex pattern to match anything inside square brackets '[...]'
    pattern = re.compile(r'\[(.*?)\]')

    # Apply the transformation using the nested reorder_attributes function
    output_newick = pattern.sub(reorder_attributes, input_newick)

    return output_newick


def print_all_trees_with_polytomies(
        trees_with_polytomies:list[TreePolytomies], distance_pairs: pandas.Series, verbose=True
) -> None:
    for idx, tp in enumerate(trees_with_polytomies):
        print(f"{idx = }")
        X: list[int] = tp.get_nodes_with_polytomies()
        print(tp)
        if verbose:
            for x in X:
                Y: list[int] = tp.get_ys(x)
                C: list[list[str]] = [tp.get_cluster(x, y_i) for y_i in Y]
                D, _, info = dms.compute_distance_matrix(distance_pairs, C, Y)
                print(info)

        print(f"{'-'*80}\n")


def extract_file_name_from_newick(newick):
    """
    Extracts, for example, the '5_10_4_2' part from a Newick string with a specific pattern.

    Args:
        newick (str): The Newick string containing patterns like 'noD_5_10_4_2|...'.

    Returns:
        str: The extracted '5_10_4_2' value, or None if not found.
    """
    # Regular expression to match the pattern 'noD_' followed by numbers separated by underscores
    pattern = r'noD_(\d+_\d+_\d+_\d+)\|'

    # Search for the pattern in the Newick string
    match = re.search(pattern, newick)
    if match:
        return match.group(1)  # Return the captured group (5_10_4_2)
    else:
        return None  # Return None if no match is found


def read_newick_from_file(base_path, file_name):
    """
    Reads a Newick string from the specified file. The folder name is deduced
    from the file name (ga_b_c_d.pruned.tree format).

    Args:
        base_path (str): Path to the root folder (e.g., 'input/true_gene_trees/').
        file_name (str): Name of the file in the format 'ga_b_c_d.pruned.tree'.

    Returns:
        str: The Newick string from the file, or None if the file is not found or invalid.
    """
    # Validate and parse the file name using regex
    pattern = re.compile(r"^g(\d+)_(\d+)_(\d+)_(\d+)\.pruned\.tree$")
    match = pattern.match(file_name)

    if not match:
        print(f"Invalid file name format: {file_name}")
        return None

    # Extract a, b, c, d from the file name
    a, b, c, d = match.groups()

    # Deduce the folder name
    folder_name = f"{a}_{b}_{c}_"

    # Construct the full file path
    file_path = os.path.join(base_path, folder_name, file_name)

    # Check if the file exists
    if not os.path.isfile(file_path):
        print(f"File not found: {file_path}")
        return None

    # Read the Newick string from the file
    with open(file_path, 'r') as f:
        real_newick = f.read().strip()

    return real_newick


# Toño's advice to work with the triplets part:
# Convert to nxTree
# -----------------
def custom_tree(nhx):
    T= read_nhxx(nhx)
    for x in T:
        T.nodes[x]['event']= 'S'
        if T.out_degree(x) == 0:  # TODO: ask Toño if it is valid; CHANGES a leaf from 'noD_5_10_4_2|G17_7' to 'G17_7'
            T.nodes[x]['label'] = T.nodes[x]['label'].split('|')[-1]
    return T


######################
# Compue performance #
######################
def triple_performance(tree: nx.DiGraph, real_tree: nx.DiGraph):
    # Todo: should the 'real_tree' be the first or second argument?
    tree1_trples= set(get_triplets(tree, color='label'))
    tree2_trples= set(get_triplets(real_tree, color='label'))
    tree2_sets= set(map(frozenset, tree2_trples))

    TP= tree1_trples . intersection( tree2_trples )
    FP= tree1_trples - tree2_trples
    FN= tree2_trples - tree1_trples
    C= {X for X in FP if frozenset(X) in tree2_sets} # Contradictory triples

    return len(TP), len(FP), len(FN), len(C)


def get_triplets(tree, event='event', color= 'color', root_event= 'S', loss_leafs= 'X'):
    """
    return a tuple (a,b,c), where a and b are the ingroup
    and c is the outgroup.
    """
    I= {} # Dictionary for induced leafs
    for x in nx.dfs_postorder_nodes(tree):
        if tree.out_degree(x) == 0:
            if tree.nodes[x][color] == loss_leafs:
                I[x]= {  }
            else:
                I[x]= { tree.nodes[x][color] }
        else:
            if x!=tree.root:
                I[x]= set( chain.from_iterable((I[x1] for x1 in tree[x])) )
            if tree.nodes[x][event] == root_event:
                for triple in _get_triples_from_root(tree, x, I):
                    yield triple


def _get_triples_from_root(tree, node, I):
    for x0, x1 in combinations(tree[node], 2):
        for triple in chain(_get_triplets_from_groups(tree, x0, x1, I),
                            _get_triplets_from_groups(tree, x1, x0, I)):
            yield triple


def _get_triplets_from_groups(tree, x_out, x_in, I):
    P= combinations( I[x_in], 2 )
    for (a,(b,c)) in product( I[x_out], P ):
        if len({a,b,c})==3:
            yield tuple(sorted((b,c)))+(a,)


def get_precision_recall_contradiction(tree: nx.DiGraph, real_tree: nx.DiGraph) -> tuple[float, float, float]:
    """
    Computes precision, recall, and contradiction based on the comparison of two trees.

    Args:
        tree (nx.DiGraph): The tree to evaluate.
        real_tree (nx.DiGraph): The real (ground-truth) tree.

    Returns:
        tuple[float, float, float]: Precision, recall, and contradiction values.
                                      NaN is returned for undefined values (e.g., division by zero).
    """
    tp, fp, fn, contradictory = triple_performance(tree, real_tree)

    # Compute precision, recall, and contradiction with safe division
    precision: float = tp / (tp + fp) if (tp + fp) > 0 else math.nan
    recall: float = tp / (tp + fn) if (tp + fn) > 0 else math.nan
    contradiction: float = contradictory / (tp + fn) if (tp + fn) > 0 else math.nan

    return precision, recall, contradiction


def is_diagonal_zero_and_nan_elsewhere(D: np.ndarray) -> bool:
    """
    Checks if a matrix has zeros on the diagonal and NaN everywhere else.

    Args:
        D (np.ndarray): The input matrix.

    Returns:
        bool: True if the matrix satisfies the condition, False otherwise.
    """

    # Check if diagonal elements are all zeros
    diagonal_zero = np.all(np.diag(D) == 0)

    # Create a mask for the off-diagonal elements
    off_diagonal_mask = ~np.eye(D.shape[0], dtype=bool)  # True for off-diagonal elements

    # Check if off-diagonal elements are all NaN
    off_diagonal_nan = np.all(np.isnan(D[off_diagonal_mask]))

    return diagonal_zero and off_diagonal_nan


def load_distance_pairs_and_trees_with_polytomies(
        hits_path: str, trees_path: str
) -> tuple[pandas.Series, list[TreePolytomies]]:
    # Load the hits and gtrees data from input files
    distance_pairs = load_hits_compute_distance_pairs(hits_path)  # Load distances
    gTrees = read_csv(trees_path, sep='\t')                             # Load trees
    gTrees = gTrees.set_index('OG').tree.apply(read_nhxx)               # Load trees
    trees_with_polytomies = get_trees_with_polytomies(gTrees)     # Identify those trees with polytomies

    return distance_pairs, trees_with_polytomies
