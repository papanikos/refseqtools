#!/usr/bin/env python

from collections import defaultdict
import argparse
from utils import *

parser = argparse.ArgumentParser(description=
                                 "Given a pre-processed RefSeq release catalog file with 3 fields (accession, taxid, size) "
                                 "create a json file that contains sequence content information for each taxid. "
                                 "This is of the form { taxid: { acc : size, ... }, ... }"
                                 
                                 "The original RefSeq catalog is available on: "
                                 "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/RefSeq-release<version>.catalog.gz"
                                 )
parser.add_argument('-i',
                    dest='catalog_in',
                    help="Path to the refseq catalog",
                    required=True
                    )
parser.add_argument('-o',
                    dest='json_out',
                    help="Path to the json file to be written",
                    required=True
                    )


def load_acc2taxid_map(map_file, json_file):
    """
    Given a refseq catalog file, write a json file aggregating
    sequence content information per taxid
    :param map_file: Refseq catalog input file
    :param json_file: JSON output file to write
    """
    d = defaultdict(list)
    entries = []
    with optionally_compressed_handle(map_file, 'r') as fmap:
        for l in fmap:
            fields = l.split('\t')
            acc = fields[0].strip()
            taxid = fields[1].strip()
            size = int(fields[2].strip())
            entries += [(acc, taxid, size)]

    for acc, taxid, size in entries:
        if taxid not in d:
            d[taxid] = {acc: size}
        else:
            d[taxid].update({acc: size})

    with open(json_file, 'w') as fj:
        json.dump(d, fj)


if __name__ == '__main__':
    args = parser.parse_args()
    load_acc2taxid_map(args.catalog_in, args.json_out)


