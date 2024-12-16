import src.Utils.Utils as utils


if __name__ == "__main__":
    # Input Newick string with mixed-order attributes
    input_newick = "(([species=H7;node_id=16;label=noD_5_10_4_2|G17_7],[node_id=18;label=X])[node_id=17;label=S],(([species=H4;node_id=11;label=noD_5_10_4_2|G13_4],([species=H2;node_id=12;label=noD_5_10_4_2|G11_2],[species=H1;node_id=10;label=noD_5_10_4_2|G10_1])[node_id=13;label=S])[node_id=14;label=S],(([species=H4;node_id=1;label=noD_5_10_4_2|G8_4],([species=H2;node_id=0;label=noD_5_10_4_2|G6_2],[species=H1;node_id=2;label=noD_5_10_4_2|G5_1])[node_id=3;label=S])[node_id=4;label=S],([species=H4;node_id=6;label=noD_5_10_4_2|G3_4],([species=H2;node_id=7;label=noD_5_10_4_2|G1_2],[species=H1;node_id=5;label=noD_5_10_4_2|G0_1])[node_id=8;label=S])[node_id=9;label=S]))[node_id=15;label=D])[node_id=19;label=S];"

    # Transform the input Newick
    output_newick = utils.transform_newick(input_newick)

    # Print the result
    print("Transformed Newick:")
    print(output_newick)
