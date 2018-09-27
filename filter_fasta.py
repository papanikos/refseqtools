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


refseq_genomic_prefixes = {"NC", "NW", "NG", "NT", "NZ_CP", "NZ_CM", "NZ", "AC"}

class FastaFilterer:
    """
    An object that filters sequences from a given fasta file, based on given taxonomic identifiers
    or accession prefixes.
    """

    def __init__(self, fp_in,
                 include_taxa, exclude_taxa,
                 include_accessions, exclude_accessions,
                 accessions_infile, accessions_exfile,
                 acc2taxid_json):

        self.fp_in = fp_in

        # Load the taxonomy info if taxa based filtering is defined
        if include_taxa or exclude_taxa:
            self.taxonomy_filtering = True
            print("Loading taxonomy information")
            self.acc2taxid = get_acc2taxid_map(acc2taxid_json)
            self.ncbiTree = NCBITaxa()
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
        if self.taxonomy_filtering and self.accession_prefix_filtering and self.accession_file_filtering:
            taxonomy_accessions = self.create_taxonomy_accessions_set()
            if self.tmode in ["both", "inclusive_only"]:
                tmp_accessions = fasta_accessions.intersection(taxonomy_accessions)
                if self.pmode in ["both", "inclusive_only"]:
                    final_accessions = filter_accession_list_on_prefix(tmp_accessions,
                                                                       self.acc_prefixes,
                                                                       exclusive = False)
                elif self.pmode == "exclusive_only":
                    final_accessions = filter_accession_list_on_prefix(tmp_accessions,
                                                                       self.acc_prefixes,
                                                                       exclusive = True)
                else:
                    print("Why am I here...?")
            elif self.tmode == "exclusive_only":
                tmp_accessions = fasta_accessions.difference(taxonomy_accessions)
                if self.pmode in ["both", "inclusive_only"]:
                    final_accessions = filter_accession_list_on_prefix(tmp_accessions,
                                                                       self.acc_prefixes,
                                                                       exclusive = False)
                elif self.pmode == "exclusive_only":
                    final_accessions = filter_accession_list_on_prefix(tmp_accessions,
                                                                       self.acc_prefixes,
                                                                       exclusive = True)
                else:
                    print("Why am I here")
        # 2. FILE AND PREFIX
        elif self.accession_file_filtering and self.accession_prefix_filtering:
            if self.amode in ["both", "inclusive_only"]:
                tmp_accessions = fasta_accessions.intersection(self.pre_accessions)
                if self.pmode in ["both", "inclusive_only"]:
                    final_accessions = filter_accession_list_on_prefix(tmp_accessions,
                                                                       self.acc_prefixes,
                                                                       exclusive=False)
                elif self.pmode == "exclusive_only":
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
                    final_accessions = filter_accession_list_on_prefix(tmp_accessions,
                                                                       self.acc_prefixes,
                                                                       exclusive=True)
            else:
                print("Why am I here...")
        # 3. FILE ONLY
        elif self.accession_file_filtering:
            if self.amode in ["both", "inclusive_only"]:
                final_accessions=fasta_accessions.union(self.pre_accessions)
            elif self.amode == "exclusive_only":
                final_accessions = fasta_accessions.difference(self.pre_accessions)

        # 4. PREFIX ONLY
        elif self.accession_prefix_filtering:
            if self.pmode in ["both", "inclusive_only"]:
                final_accessions = filter_accession_list_on_prefix(fasta_accessions,
                                                                   self.acc_prefixes,
                                                                   exclusive=False)
            elif self.pmode == "exclusive_only":
                final_accessions = filter_accession_list_on_prefix(fasta_accessions,
                                                                   self.acc_prefixes,
                                                                   exclusive=True)
        # 5. TAXONOMY AND PREFIX
        elif self.taxonomy_filtering and self.accession_prefix_filtering:
            taxonomy_accessions = self.create_taxonomy_accessions_set()
            if self.tmode in ["both", "inclusive_only"]:
                tmp_accessions = fasta_accessions.intersection(taxonomy_accessions)
                if self.pmode in ["both", "inclusive_only"]:
                    final_accessions = filter_accession_list_on_prefix(tmp_accessions,
                                                                       self.acc_prefixes,
                                                                       exclusive = False)
                elif self.pmode == "exclusive_only":
                    final_accessions = filter_accession_list_on_prefix(tmp_accessions,
                                                                       self.acc_prefixes,
                                                                       exclusive = True)

        # 6. TAXONOMY AND ACCESSION FILE
        elif self.taxonomy_filtering and self.accession_file_filtering:
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

    def make_output_filename(self, fp_in, fasta_type="genomic"):
        filtered_prefix = "filtered_" + fasta_type
        self.fp_out = fp_in.replace(fasta_type, filtered_prefix)
        return self.fp_out

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
                    counter += 1
        else:
            print("No sequences were matched")

        return counter, fp_out

    ### THIS IS OBSOLETE
    def filter_accession_to_file(self, fp_in, acc_prefixes):
        fp_out = self.make_output_filename(fp_in)
        counter = 0
        with optionally_compressed_handle(fp_out, 'w') as fout, optionally_compressed_handle(fp_in, 'r') as fin:
            input_seq_iterator = SeqIO.parse(fin, format="fasta")
            for record in input_seq_iterator:
                if not self.pre_accessions:
                    if keep_accession(record.id, acc_prefixes):
                        counter += 1
                        SeqIO.write(record, fout, format="fasta")
                elif keep_accession(record.id, self.acc_prefixes) or record.id in self.pre_accessions:
                    counter += 1
                    SeqIO.write(record, fout, format="fasta")
        if counter == 0:
            print("No accessions were found matching the parameters.")
            os.remove(fp_out)
        return counter

    def write_accessions_to_file(self, fp_in, accessions):
        fp_out = self.make_output_filename(fp_in)
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
        acc_ex = [acc.strip() for acc in args.acc_ex.split(',') if args.acc_ex]
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

    filterObj = FastaFilterer(fp_in=args.fasta_in,
                              acc2taxid_json=args.json_fp,
                              include_taxa=taxa_in, exclude_taxa=taxa_ex,
                              include_accessions=acc_in, exclude_accessions=acc_ex,
                              accessions_infile=acc_in_file, accessions_exfile=acc_ex_file
                             )

    # Create the set of fasta accessions
    fasta_accessions = filterObj.get_fasta_accessions_set()
    # Filter the accessions based on input criteria
    final_accessions = filterObj.create_final_accessions_set()
    # if filtering didn't change the original set
    if final_accessions == fasta_accessions:
        # For now just rename the file - mark as processed
        outfile = filterObj.make_output_filename(filterObj.fp_in)
        os.rename(filterObj.fp_in, outfile)
        print("{}\t(renamed)\t{}\t{}\t{}".format(filterObj.fp_in,
                                                 len(final_accessions),
                                                 outfile,
                                                 len(final_accessions)))
    elif len(final_accessions) == 0:
        os.remove(filterObj.fp_in)
        print("{}\t(removed)\t{}\tNA\t0".format(filterObj.fp_in, len(fasta_accessions)))
    else:
        fp_in, seqs_in, fp_out, seqs_out = filterObj.write_accessions_to_file(filterObj.fp_in, final_accessions)
        os.remove(filterObj.fp_in)
        print("{}\t(removed)\t{}\t{}\t{}".format(fp_in, seqs_in, fp_out, seqs_out))

    print("Done!")

