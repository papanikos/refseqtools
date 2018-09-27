#!/usr/bin/env python

from ete3 import NCBITaxa
from Bio import SeqIO
import utils as U
import argparse
from annotate_tree import *
import sqlite3

parser = argparse.ArgumentParser(description=
                                 "Plot a tree with sequence content information from the Refseq releas for a given "
                                 "taxid."
                               )

parser.add_argument('-i',
                    dest='fasta_in',
                    help="Input file to filter",
                    required=True
                    )
parser.add_argument('-acc_db',
                    dest='acc_db',
                    help='sqlite3 db file path',
                    required=True
                    )
parser.add_argument('-prune',
                    dest='prune_level',
                    type=str,
                    help='Level on which the tree will be pruned. (e.g. family)'
                         'Make sure this is at least a level deeper than the given taxid',
                    required=False)
parser.add_argument('-r',
                    dest='render',
                    help='Path to file the image is to be written. Valid extensions are .svg, .pdf, .png '
                         'and are detected automatically',
                    required=False
                    )


def get_accession_list_from_fasta(fasta_file):
    accessions = []
    with U.optionally_compressed_handle(fasta_file, 'r') as fin:
        for record in SeqIO.parse(fin, format="fasta"):
            accessions.append(record.id)
    return accessions

def make_accessions_dict(accessions_list, seq_map):
    new_dict = {}
    for acc in accessions_list:
        for taxid in seq_map:
            if acc in seq_map[taxid].keys():
                if taxid not in new_dict:
                    new_dict[taxid] = {acc : seq_map[taxid][acc]}
                elif acc not in new_dict[taxid].keys():
                    new_dict[taxid][acc] = seq_map[taxid][acc]

    return new_dict

def get_accession_taxids_from_sqldb(acc_list, sql_db):
    result_dic = {}
    conn = sqlite3.connect(sql_db)
    c = conn.cursor()
    placeholder = '?'
    acc_batches = [acc_list[i:i+500] for i in range(0, len(acc_list), 500)]
    print("split {} accessions in {} batches of 500".format(len(acc_list), len(acc_batches)))
    batch_counter = 0
    for batch in acc_batches:
        batch_counter+=1
        print("processing batch {}".format(batch_counter))
        placeholders = ', '.join(placeholder for _ in batch)
        query = "SELECT acc,taxid,size FROM acc_sizes WHERE acc in (%s)" % placeholders
        c.execute(query, batch)
        result = c.fetchall()
        for acc, taxid, size in result:
            if taxid not in result_dic:
                result_dic[taxid] = {acc: size}
            else:
                result_dic[taxid].update({acc: size})
    conn.close()

    return result_dic


if __name__ == '__main__':
    args = parser.parse_args()
    print("Getting accession IDs")
    acc_list = get_accession_list_from_fasta(args.fasta_in)
    print("Reducing sequences to selected accessions")
    filtered_dict = get_accession_taxids_from_sqldb(acc_list, args.acc_db)
    ncbi = NCBITaxa()
    print("Calculateing tree")
    tree = ncbi.get_topology(list(filtered_dict.keys()))
    for n in tree.traverse():
        n.add_features(seq_content = 0)
    for n in tree.traverse():
        n.seq_content = sum_node_sequences(n, filtered_dict, keep_accessions=None)

    if args.prune_level:
        print("Pruning tree")
        final_nodes = nodes_to_prune(tree, args.prune_level)
        tree.prune(final_nodes)
    print("Plotting")
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.layout_fn = layout
    if args.render:
        tree.render(args.render, dpi=300, tree_style=ts)
    else:
        tree.show(tree_style=ts)
