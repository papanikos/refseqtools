#!/usr/bin/env python

import os
import gzip
from Bio import SeqIO
import json
import ftplib
import itertools
import re


def check_or_create_dir(dirname):
    """
    Check if a directory exists. Create it if it doesn't
    :param dirname: Directory name
    :return:
    """
    if not os.path.isdir(dirname):
        os.mkdir(dirname)

def single_column_file_to_list(filename):
    data_list = []
    with open(filename, 'r') as fin:
        if len(fin.readline().split()) > 1:
            raise IOError("Input file {} should contain values in a single column".format(filename))
        else:
            fin.seek(0)
            for l in fin:
                data_list.append(l.strip())
    return data_list


def select_genomic_fnas(files_list):
    """
    Given a list of files select only the ones that contain
    end with genomic.fna.gz
    :param files_list:
    :return: list of files
    """
    fna_files = []
    for f in files_list:
        if f.endswith(".genomic.fna.gz"):
            fna_files += [f]
    return fna_files


def is_gz(path):
    """
    Return true if gzipped file
    :param path: path to file
    :return: boolean
    """
    return path.endswith(".gz") or path.endswith(".z")


def optionally_compressed_handle(path, mode):
    """
    Return a file handle that is optionally gzip compressed
    :param path: path
    :param mode: mode
    :return: handle
    """
    if mode == "r" or mode == "rb":
        mode = "rt"
    if mode == "w" or mode == "wb":
        mode = "wt"
    if is_gz(path):
        return gzip.open(path, mode=mode)
    else:
        return open(path, mode=mode)


def get_acc2taxid_map(json_fp):
    """
    Read a json file into a dictionary
    :param json_fp: Path to json file
    :return:
    """
    with open(json_fp, 'r') as fin:
        mapping = json.load(fin)
    return mapping


def create_full_taxa_list(taxids, ncbi_tree, include_parent=True):
    """
    Get all the descendants for a list of taxonomy ids.
    :param taxids: A list of taxonomy ids to get the descendants for
    :param ncbi_tree: An ete3.NCBITaxa() object
    :param include_parent: If true, the query taxid is included in the resulting list
    :return: A list taxids
    """

    taxa = []
    for taxid in taxids:
        descendants = ncbi_tree.get_descendant_taxa(taxid, intermediate_nodes=True)
        taxa += [descendants]

    if include_parent:
        taxa.append(taxids)

    if len(taxa) > 1:
        taxa_flat = [taxon for sublist in taxa for taxon in sublist]
        return taxa_flat
    else:
        return taxa


def filter_accession_list_on_prefix(acc_list, acc_prefixes, exclusive = False):
    """
    Given a list of accessions and a list of prefixes, keep the accessions
    starting with `prefix`. For NZ accessions NZ_CP and NZ_CM suffixes are checked separately.
    :param acc_list: A list of RefSeq accessions, i.e. containing '_'
    :param acc_prefixes: A list of prefixes. Should be in ['AC', 'NC', 'NG', 'NT', 'NW', 'NZ', 'NZ_CP', 'NZ_CM']
    :return:
    """
    filtered_accessions = []
    for acc in acc_list:
        if keep_accession(acc, acc_prefixes) and not exclusive:
            filtered_accessions.append(acc)
        elif not keep_accession(acc,acc_prefixes) and exclusive:
            filtered_accessions.append(acc)
    return filtered_accessions

def keep_accession(acc, acc_prefixes):
    fields = acc.split('_')
    CP_reference = re.compile('CP(?=[0-9])')
    CM_reference = re.compile('CM(?=[0-9])')
    if fields[0] == 'NZ':
        if 'NZ_CP' in acc_prefixes and re.match(CP_reference, fields[1]) is not None:
            return True
        elif 'NZ_CM' in acc_prefixes and re.match(CM_reference, fields[1]) is not None:
            return True
        elif 'NZ' in acc_prefixes:
            return True
    elif fields[0] in acc_prefixes:
        return True

class FastaFilterer:
    """
    An object that filters sequences from a given fasta file, based on given taxonomic identifiers.
    """
    def __init__(self, fp_in,
                 include_taxa, exclude_taxa,
                 include_accessions, exclude_accessions,
                 acc2taxid, ncbiTree):

        self.fp_in = fp_in
        self.acc2taxid = acc2taxid
        self.ncbiTree = ncbiTree

        if include_taxa and exclude_taxa:
            taxa_in = set(create_full_taxa_list(include_taxa, self.ncbiTree, include_parent=True))
            taxa_ex = set(create_full_taxa_list(exclude_taxa, self.ncbiTree, include_parent=True))
            if taxa_ex.issuperset(taxa_in):
                self.final_taxa_list = list(taxa_in)
            else:
                self.final_taxa_list = list(taxa_in - taxa_ex)

        else:
            self.final_taxa_list = create_full_taxa_list(include_taxa, self.ncbiTree, include_parent=True)

        if include_accessions and exclude_accessions:
            self.acc_prefixes = list(set(include_accessions) - set(exclude_accessions))
        else:
            self.acc_prefixes = include_accessions

    def create_list_of_accessions(self):
        accessions = []
        for taxid in self.final_taxa_list:
            tmp_ac = self.acc2taxid.get(str(taxid), None)
            if tmp_ac is not None:
                accessions += [list(tmp_ac.keys())]
        accessions_flat = [accession for sublist in accessions for accession in sublist]
        final_accessions = filter_accession_list_on_prefix(accessions_flat, self.acc_prefixes)
        return final_accessions

    def make_output_filename(self, fp_in, fasta_type="genomic"):
        filtered_prefix = "filtered_" + fasta_type
        fp_out = fp_in.replace(fasta_type, filtered_prefix)
        return fp_out

    def peek(self, iterable):
        """
        This checks if a generator is empty
        taken from https://stackoverflow.com/a/664239
        :param iterable: A generator
        :return:
        """
        try:
            first = next(iterable)
        except StopIteration:
            return None
        return first, itertools.chain([first], iterable)

    def filter_records_to_file(self, fp_in, accessions_list):
        counter = 0
        fp_out = self.make_output_filename(fp_in)
        # Create an iterator
        input_seq_iterator = SeqIO.parse(optionally_compressed_handle(fp_in, 'r'), "fasta")
        # Filter the sequences based on the ids
        keepers = (record for record in input_seq_iterator
                   if record.id in accessions_list)

        # Check if the iterator is empty so empty files won't be open for writing
        res = self.peek(keepers)
        if res is not None:
            first, seqs = res
            with optionally_compressed_handle(fp_out, 'w') as fout:
                for record in seqs:
                    SeqIO.write(record, fout, format="fasta")
        else:
            print("No sequences were matched")

        return counter, fp_out

    
class NcbiFTPConnector:
    def __init__(self, ftp_address="ftp.ncbi.nlm.nih.gov",
                 release_dir= "refseq/release"):
        
        # print("Connecting to {}".format(ftp_address))
        
        ftp = ftplib.FTP(ftp_address)
        
        try:
            ftp.login()
        except Exception as e:
            raise e('Connection failed')

        self.ftp = ftp
        self.ftp_address = ftp_address
        self.release_dir = release_dir
        
    def go_to_dir(self, target_dir):
        target_dir = os.path.join(self.release_dir, target_dir)
        self.ftp.cwd(target_dir)
