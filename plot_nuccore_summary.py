#!/usr/bin/env python3

from ete3 import NCBITaxa, faces, TreeStyle
import re
import sys
from annotate_tree import *


# Each downloaded entry is numbered
entry = "^[0-9]*\."
pattern = re.compile(entry)


def parse_nuccore_summary(summary_txt):
    d = {}
    with open(summary_txt, 'r') as fin:
        for l in fin:
            if re.match(pattern, l):
                fields = l.split(' ')
                species = ' '.join(fields[1:3]).rstrip()
                ## Handle specific errors that occur
                if species == "Drosophila busckii,":
                    species = "Drosophila busckii"
                if species == "Malus x":
                    species = "Malus x domestica"
                # Move to the next line that contains the sequence info
                next_line = next(fin)
                size = int(next_line.split(' ')[0].replace(',', ''))

                if species not in d:
                    d[species] = size
                else:
                    d[species] += size
    return d


def merge_dics(nuccore_dic, taxids_dic):
    """
    Merge a dictionary that contains sequence sizes for species and
    a dictionary with species - taxid mapping.
    :param nuccore_dic: A dictionary that contains the total sequence size for a species,
    after parsing the nuccore summary.txt file
    :param taxids_dic: A dictionary that contains the mapping of each human readable species name to its taxid
    :return:
    """
    new_dic = {}
    for species in nuccore_dic.keys():
        new_dic[taxids_dic[species][0]] = {species: nuccore_dic[species]}
    return new_dic


if __name__ == '__main__':
    # {species : sequence_size }
    nuccore = parse_nuccore_summary(sys.argv[1])
    species_list = list(nuccore.keys())

    ncbi = NCBITaxa()

    # {species : taxid }
    taxids_dic = ncbi.get_name_translator(species_list)

    taxids_list = [i[0] for i in list(taxids_dic.values())]

    tree = ncbi.get_topology(taxids_list)

    all_info = merge_dics(nuccore, taxids_dic)

    for node in tree.traverse():
        node.add_features(seq_content=0)

    for node in tree.traverse():
        sum_node_sequences(node, all_info)

    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.layout_fn = layout
    tree.show(tree_style=ts)

    dir, fname = os.path.split(sys.argv[1])
    img_file_fp = os.path.join(dir, fname.replace("txt", "png"))
    tree.render(img_file_fp, dpi=300, tree_style=ts)

