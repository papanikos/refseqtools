#!/usr/bin/env python

from Bio import SeqIO
from utils import optionally_compressed_handle
import sys

def fasta_to_dict(fasta_in):
    with optionally_compressed_handle(fasta_in, 'r') as fin:
        seq_dict = SeqIO.to_dict(SeqIO.parse(fin, "fasta"))
    print(seq_dict)
    return seq_dict

def write_sorted_seq_dict(seq_dict, fasta_out):
    with optionally_compressed_handle(fasta_out, 'w') as fout:
        sorted_keys = sorted(seq_dict.keys())
        for id in sorted_keys:
            SeqIO.write(seq_dict[id], fout, format='fasta')

if __name__ == '__main__':
    dic = fasta_to_dict(sys.argv[1])
    write_sorted_seq_dict(dic, sys.argv[2])

