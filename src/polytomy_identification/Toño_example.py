from pandas import read_csv
from revolutionhtl.nhxx_tools import read_nhxx, get_nhx
from revolutionhtl.nxTree import induced_colors

# Parameters
trees_path = '../../input/tl_project.reconciliation.tsv'

# Functions


def is_polytomi(T, node):
    return len(T[node]) > 2


def get_polytomies(T):
    return [node for node in T if is_polytomi(T, node)]


if __name__ == "__main__":
    # Load trees
    gTrees = read_csv(trees_path, sep='\t')
    gTrees = gTrees.set_index('OG').tree.apply(read_nhxx)

    # Identify polytomies
    polytomies = gTrees.apply(get_polytomies)

    # Example of how to get cluster of leaves
    example_tree = gTrees[0]
    example_node = 1
    induced_colors(example_tree, example_node, 'label')

    # Example of how to convert back to newick
    nhx = get_nhx(example_tree, name_attr='label')
    print(nhx)
