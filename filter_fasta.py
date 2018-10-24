from ete3 import NCBITaxa
import argparse
from utils import *

parser = argparse.ArgumentParser(description=
                                 "Filter fasta files based on taxonomy and accession types"
                                 )
parser.add_argument('-i',
                    dest='fasta_in',
                    help="Input file to filter",
                    required=True
                    )
parser.add_argument('-o',
                    dest="fasta_out",
                    help="Fasta ouput file. Depending on extensions can be gzipped or not",
                    required=True)
parser.add_argument('--include-taxa',
                    dest='taxa_in',
                    help="A comma separated list of taxa to be filtered",
                    required=False
                    )
parser.add_argument('--exclude-taxa',
                    dest='taxa_ex',
                    help='Taxa to be excluded',
                    required=False
                    )
parser.add_argument('--include-accessions',
                    dest='acc_in',
                    help='List of accession types to be included',
                    required=False
                    )
parser.add_argument('--exclude-accessions',
                    dest='acc_ex',
                    help='List of accession types to be excluded',
                    required=False
                    )
parser.add_argument('-j',
                    dest='json_fp',
                    help='JSON file with the mapping of taxonomy ids to accessions',
                    required=True
                    )
parser.add_argument('--include-accessions-from-file',
                    dest='acc_in_file',
                    help="A file containing accessions to be included, one per line",
                    required=False)
parser.add_argument('--exclude-accessions-from-file',
                    dest='acc_ex_file',
                    help="A file containing accessions to be excluded, one per line",
                    required=False)
parser.add_argument("--dbfile",
                    dest="dbfile",
                    help="Path to instantiated ete3 sqlite db",
                    required=False)


refseq_genomic_prefixes = {"NC", "NW", "NG", "NT", "NZ_CP", "NZ_CM", "NZ", "AC"}


class FastaFilterer:
    """
    An object that filters sequences from a given fasta file, based on given taxonomic identifiers
    or accession prefixes.
    """

    def __init__(self, fp_in, fp_out,
                 include_taxa, exclude_taxa,
                 include_accessions, exclude_accessions,
                 accessions_infile, accessions_exfile,
                 acc2taxid_json, dbfile):

        self.fp_in = fp_in
        self.fp_out = fp_out

        # Load the taxonomy info if taxa based filtering is defined
        if include_taxa or exclude_taxa:
            self.taxonomy_filtering = True
            print("Loading taxonomy information")
            self.acc2taxid = get_acc2taxid_map(acc2taxid_json)
            self.ncbiTree = NCBITaxa(dbfile=dbfile)
            if include_taxa and exclude_taxa:
                self.tmode = "both"
                taxa_in = set(create_full_taxa_list(include_taxa, self.ncbiTree, include_parent=True))
                taxa_ex = set(create_full_taxa_list(exclude_taxa, self.ncbiTree, include_parent=True))
                if taxa_ex.issuperset(taxa_in):
                    self.final_taxa_list = list(taxa_in)
                else:
                    self.final_taxa_list = list(taxa_in - taxa_ex)

            elif include_taxa and not exclude_taxa:
                self.tmode = "inclusive_only"
                self.final_taxa_list = create_full_taxa_list(include_taxa, self.ncbiTree, include_parent=True)
            elif exclude_taxa and not include_taxa:
                self.tmode = "exclusive_only"
                self.final_taxa_list = create_full_taxa_list(exclude_taxa, self.ncbiTree, include_parent=True)
        else:
            self.final_taxa_list = None
            self.taxonomy_filtering = None

        if include_accessions or exclude_accessions:
            self.accession_prefix_filtering = True
            if include_accessions and exclude_accessions:
                self.pmode = "both"
                self.acc_prefixes = list(set(include_accessions) - set(exclude_accessions))
            elif include_accessions and not exclude_accessions:
                self.pmode = "inclusive_only"
                self.acc_prefixes = include_accessions
            elif exclude_accessions and not include_accessions:
                self.pmode = "exclusive_only"
                self.acc_prefixes = refseq_genomic_prefixes - set(exclude_accessions)
        else:
            self.acc_prefixes = None
            self.accession_prefix_filtering = None

        if accessions_infile or accessions_exfile:
            self.accession_file_filtering = True
            if accessions_infile and accessions_exfile:
                self.amode = "both"
                faccessions_in = set(single_column_file_to_list(accessions_infile))
                faccessions_ex = set(single_column_file_to_list(accessions_exfile))
                self.pre_accessions = faccessions_in - faccessions_ex
            elif accessions_infile and not accessions_exfile:
                self.amode = "inclusive_only"
                self.pre_accessions = set(single_column_file_to_list(accessions_infile))
            elif accessions_exfile and not accessions_infile:
                self.amode = "exclusive_only"
                self.pre_accessions = set(single_column_file_to_list(accessions_exfile))
        else:
            self.pre_accessions = None
            self.accession_file_filtering = None

    def create_taxonomy_accessions_set(self):
        accessions = []
        for taxid in self.final_taxa_list:
            tmp_ac = self.acc2taxid.get(str(taxid), None)
            if tmp_ac is not None:
                accessions += [list(tmp_ac.keys())]
        accessions_flat = [accession for sublist in accessions for accession in sublist]
        taxonomy_accessions = filter_accession_list_on_prefix(accessions_flat, self.acc_prefixes)
        return set(taxonomy_accessions)

    def get_fasta_accessions_set(self):
        seq_iterator = SeqIO.parse(optionally_compressed_handle(self.fp_in, 'r'), format="fasta")
        fasta_accessions = {record.id for record in seq_iterator}
        return fasta_accessions

    def create_final_accessions_set(self):
        fasta_accessions = self.get_fasta_accessions_set()
        # 1. ALL FILTERS ENABLED
        # FIX ME
        # The self.amode is not covered here...
        if self.taxonomy_filtering \
                and self.accession_prefix_filtering \
                and self.accession_file_filtering:

            taxonomy_accessions = self.create_taxonomy_accessions_set()
            if self.tmode in ["both", "inclusive_only"]:
                tmp_accessions = fasta_accessions.intersection(taxonomy_accessions)

                if self.pmode in ["both", "inclusive_only"]:
                    final_accessions = filter_accession_list_on_prefix(tmp_accessions,
                                                                       self.acc_prefixes,
                                                                       exclusive=False)
                elif self.pmode == "exclusive_only":
                    final_accessions = filter_accession_list_on_prefix(tmp_accessions,
                                                                       self.acc_prefixes,
                                                                       exclusive=True)
                else:
                    print("Why am I here...?")

            elif self.tmode == "exclusive_only":
                tmp_accessions = fasta_accessions.difference(taxonomy_accessions)
                if self.pmode in ["both", "inclusive_only"]:
                    final_accessions = filter_accession_list_on_prefix(tmp_accessions,
                                                                       self.acc_prefixes,
                                                                       exclusive=False)
                elif self.pmode == "exclusive_only":
                    print("excluding prefixes")
                    final_accessions = filter_accession_list_on_prefix(tmp_accessions,
                                                                       self.acc_prefixes,
                                                                       exclusive=True)
                else:
                    print("Why am I here")

        # 2. FILE AND PREFIX
        elif self.accession_file_filtering \
                and self.accession_prefix_filtering \
                and not self.taxonomy_filtering:

            if self.amode in ["both", "inclusive_only"]:
                tmp_accessions = fasta_accessions.intersection(self.pre_accessions)
                if self.pmode in ["both", "inclusive_only"]:
                    print("including prefixes")
                    final_accessions = filter_accession_list_on_prefix(tmp_accessions,
                                                                       self.acc_prefixes,
                                                                       exclusive=False)
                elif self.pmode == "exclusive_only":
                    print("excluding prefixes")
                    final_accessions = filter_accession_list_on_prefix(tmp_accessions,
                                                                       self.acc_prefixes,
                                                                       exclusive=True)
                else:
                    print("Why am I here")
            elif self.amode == "exclusive_only":
                tmp_accessions = fasta_accessions.difference(self.pre_accessions)
                if self.pmode in ["both", "inclusive_only"]:
                    final_accessions = filter_accession_list_on_prefix(tmp_accessions,
                                                                       self.acc_prefixes,
                                                                       exclusive=False)
                elif self.pmode == "exclusive_only":
                    print("excluding prefixes")
                    final_accessions = filter_accession_list_on_prefix(tmp_accessions,
                                                                       self.acc_prefixes,
                                                                       exclusive=True)
            else:
                print("Why am I here...")

        # 3. FILE ONLY
        elif self.accession_file_filtering \
                and not self.accession_prefix_filtering \
                and not self.taxonomy_filtering:

            if self.amode in ["both", "inclusive_only"]:
                final_accessions = fasta_accessions.union(self.pre_accessions)
            elif self.amode == "exclusive_only":
                final_accessions = fasta_accessions.difference(self.pre_accessions)

        # 4. PREFIX ONLY
        elif self.accession_prefix_filtering \
                and not self.accession_file_filtering \
                and not self.taxonomy_filtering:

            if self.pmode in ["both", "inclusive_only"]:
                final_accessions = filter_accession_list_on_prefix(fasta_accessions,
                                                                   self.acc_prefixes,
                                                                   exclusive=False)
            elif self.pmode == "exclusive_only":
                final_accessions = filter_accession_list_on_prefix(fasta_accessions,
                                                                   self.acc_prefixes,
                                                                   exclusive=True)

        # 5. TAXONOMY AND PREFIX
        elif self.taxonomy_filtering \
                and self.accession_prefix_filtering \
                and not self.accession_file_filtering:

            taxonomy_accessions = self.create_taxonomy_accessions_set()
            if self.tmode in ["both", "inclusive_only"]:
                tmp_accessions = fasta_accessions.intersection(taxonomy_accessions)
            elif self.tmode == "exclusive_only":
                tmp_accessions = fasta_accessions.difference(taxonomy_accessions)
            else:
                tmp_accessions = set()

            if self.pmode in ["both", "inclusive_only"]:
                final_accessions = filter_accession_list_on_prefix(tmp_accessions,
                                                                   self.acc_prefixes,
                                                                   exclusive=False)
            elif self.pmode == "exclusive_only":
                final_accessions = filter_accession_list_on_prefix(tmp_accessions,
                                                                   self.acc_prefixes,
                                                                   exclusive=True)

        # 6. TAXONOMY AND ACCESSION FILE
        elif self.taxonomy_filtering \
                and self.accession_file_filtering \
                and not self.accession_prefix_filtering:

            taxonomy_accessions = self.create_taxonomy_accessions_set()
            if self.tmode in ["both", "inclusive_only"]:
                tmp_accessions = fasta_accessions.intersection(taxonomy_accessions)
                if self.amode in ["both", "inclusive_only"]:
                    final_accessions = tmp_accessions.intersection(self.pre_accessions)
                elif self.amode == "exclusive_only":
                    final_accessions = tmp_accessions.difference(self.pre_accessions)

        else:
            print("Case not covered")

        return set(final_accessions)

    def write_accessions_to_file(self, fp_in, fp_out, accessions):
        out_counter, in_counter = 0, 0
        with optionally_compressed_handle(fp_out, 'w') as fout, optionally_compressed_handle(fp_in, 'r') as fin:
            input_seq_iterator = SeqIO.parse(fin, format="fasta")
            for record in input_seq_iterator:
                in_counter += 1
                if record.id in accessions:
                    out_counter += 1
                    SeqIO.write(record, fout, format="fasta")
        return fp_in, in_counter, fp_out, out_counter


if __name__ == '__main__':
    args = parser.parse_args()

    if args.taxa_in:
        taxa_in = [taxid.strip() for taxid in args.taxa_in.split(',')]
    else:
        taxa_in = None

    if args.taxa_ex:
        taxa_ex = [taxid.strip() for taxid in args.taxa_ex.split(',')]
    else:
        taxa_ex = None

    if args.acc_in:
        acc_in = [acc.strip() for acc in args.acc_in.split(',')]
    else:
        acc_in = ['AC', 'NG', 'NC', 'NW', 'NZ', 'NZ_CP', 'NZ_CM', 'NT']

    if args.acc_ex:
        acc_ex = [acc.strip() for acc in args.acc_ex.split(',')]
    else:
        acc_ex = None

    if args.acc_in_file:
        acc_in_file = args.acc_in_file
    else:
        acc_in_file = None

    if args.acc_ex_file:
        acc_ex_file = args.acc_ex_file
    else:
        acc_ex_file = None

    filterObj = FastaFilterer(fp_in=args.fasta_in, fp_out=args.fasta_out,
                              acc2taxid_json=args.json_fp,
                              include_taxa=taxa_in, exclude_taxa=taxa_ex,
                              include_accessions=acc_in, exclude_accessions=acc_ex,
                              accessions_infile=acc_in_file, accessions_exfile=acc_ex_file,
                              dbfile=args.dbfile
                             )

    # Create the set of fasta accessions
    fasta_accessions = filterObj.get_fasta_accessions_set()
    # Filter the accessions based on input criteria
    final_accessions = filterObj.create_final_accessions_set()
    # if filtering didn't change the original set
    if final_accessions == fasta_accessions:
        # For now just rename the file - mark as processed
        fp_in, seqs_in, fp_out, seqs_out = filterObj.write_accessions_to_file(filterObj.fp_in,
                                                                              filterObj.fp_out,
                                                                              final_accessions)
        print("{}\t(copied)\t{}\t{}\t{}".format(filterObj.fp_in,
                                                 len(fasta_accessions),
                                                 filterObj.fp_out,
                                                 len(final_accessions)))
    elif len(final_accessions) == 0:
        print("{}\t(untouched)\t{}\tNA\t0".format(filterObj.fp_in, len(fasta_accessions)))
    else:
        fp_in, seqs_in, fp_out, seqs_out = filterObj.write_accessions_to_file(filterObj.fp_in,
                                                                              filterObj.fp_out,
                                                                              final_accessions)
        print("{}\t(new)\t{}\t{}\t{}".format(fp_in, seqs_in, fp_out, seqs_out))

    print("Done!")
