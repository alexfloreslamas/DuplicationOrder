from pandas import read_csv
from revolutionhtl.nhxx_tools import read_nhxx
import src.Utils.Utils as utils

# Parameters
trees_path = '../input/tl_project.reconciliation.tsv'


if __name__ == "__main__":
    # Load trees
    gTrees = read_csv(trees_path, sep='\t')
    gTrees = gTrees.set_index('OG').tree.apply(read_nhxx)
    print(f"{len(gTrees) = }")

    trees_with_polytomies = utils.get_trees_with_polytomies(gTrees)
    print(f"{len(trees_with_polytomies) = }")

    print("Example:")
    tp = trees_with_polytomies[0]
    print(tp)

    # for tree in trees_with_polytomies:
    #     print(tree)
