#!/usr/bin/env python

from utils import *
import argparse


parser = argparse.ArgumentParser(description="Compare refseq releases")

parser.add_argument('-current',
                        dest='current',
                        help="current catalog",
                        required=True
                        )

parser.add_argument('-catalog',
                    dest='catalog',
                    help='Path to prevous catalog',
                    required=True
                    )

def release_to_dic(dna_catalog):
    dic={}
    with optionally_compressed_handle(dna_catalog, 'r') as fin:
        for record in fin:
            fields = record.split('\t')
            accession = fields[0].strip()
            accession_info = accession.split('.')
            accession_id, accession_version = accession_info[0], int(accession_info[1])
            taxid = fields[1].strip()
            size = fields[2]
            if accession_id not in dic:
                dic[accession_id] = accession_version
    return dic

if __name__ == '__main__':
    args = parser.parse_args()
    current_dic = release_to_dic(args.current)
    total_current = len(list(current_dic.keys()))
    print("{} in {}".format(total_current, args.current))
    previous_dic = release_to_dic(args.catalog)
    total_previous = len(list(previous_dic.keys()))
    identical = 0
    print("{} in {}".format(total_previous, args.catalog))
    different_version = 0
    only_in_current = 0
    only_in_previous = 0
    for acc in current_dic:
        if acc in previous_dic:
            if current_dic[acc] == previous_dic[acc]:
                identical += 1
            else:
                print('{}\t{}\t{}'.format(acc, current_dic[acc], previous_dic[acc]))
                different_version += 1
        else:
            only_in_current += 1


    with open("only_in_previous.txt", 'w') as fout:

        for acc in previous_dic:
            if acc not in current_dic:
                only_in_previous += 1
                fout.write("{}\n".format('.'.join([acc, str(previous_dic[acc])])))

    print("Identical: {}".format(identical))
    print("Different version: {}".format(different_version))
    print("Only in current: {}".format(only_in_current))
    print("Only in previous: {}".format(only_in_previous))

