import Utils.Utils as utils


if __name__ == "__main__":
    # Path to the hits file.
    hits_path = '../input/tl_project_alignment_all_vs_all/'
    distance_pairs = utils.load_hits_compute_distance_pairs(hits_path)

    # Example leaf identifiers
    l1: str = 'noD_5_3_1_1|G4_4'
    l2: str = 'noD_5_3_1_1|G1_1'
    l3: str = 'I do not exist'

    # Query distances
    print(utils.get_pair_distance(distance_pairs, l1, l2))    # Expected: float, say abc.def...
    print(utils.get_pair_distance(distance_pairs, l2, l1))    # Same as above, i.e, abc.def...
    print(utils.get_pair_distance(distance_pairs, l1, l3))    # Expected: np.nan
    print(utils.get_pair_distance(distance_pairs, l3, l1))    # Expected: np.nan
    print(utils.get_pair_distance(distance_pairs, l2, l3))    # Expected: np.nan
    print(utils.get_pair_distance(distance_pairs, l3, l2))    # Expected: np.nan
