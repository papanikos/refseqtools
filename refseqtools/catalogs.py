from .utils import NcbiFTPConnector, optionally_compressed_handle
from os import path


def make_output_filename(output_dir, filename):
    output_filename = path.join(output_dir, filename)
    return output_filename


def file_exists(filename):
    if path.isfile(filename):
        print("File {} already exists. Continuing...".format(filename))
        return True
    else:
        return False


def download_from_ftp(ftp_conn, filename, output_file):
    # ftp connection is already open and it doesn't close here
    if not file_exists(output_file):
        with open(output_file, 'w') as fout:
            ftp_conn.ftp.retrbinary('RETR %s' % filename, fout.write)


def get_current_release_version(output_dir, release_number_file='current_release.txt'):
    release_number_file = make_output_filename(output_dir, release_number_file)
    if not file_exists(release_number_file):
        conn = NcbiFTPConnector()
        conn.ftp.cwd('refseq/release')
        download_from_ftp(conn, 'RELEASE_NUMBER', release_number_file)
        conn.ftp.close()

    else:
        with open(release_number_file, 'r') as fin:
            current_version = int(fin.readline().strip())

    return current_version


def get_refseq_release_catalog(output_dir, release_version):
    basename = 'RefSeq-release{}.catalog.gz'.format(str(release_version))
    current_version = get_current_release_version(output_dir)
    catalog_file = make_output_filename(output_dir, basename)
    if not file_exists(catalog_file):
        conn = NcbiFTPConnector()
        if current_version == release_version:
            conn.go_to_dir('release-catalog')  # This should be stable, hence hardcoded
            print("Downloading file {}".format(basename))
            download_from_ftp(conn, basename, catalog_file)
        else:
            print("Downloading file {}".format(basename))
            conn.go_to_dir('release-catalog/archive')  # This should be stable, hence hardcoded
            download_from_ftp(conn, basename, catalog_file)

        conn.ftp.close()

    return catalog_file


def genomic_records_to_dic(catalog_file):
    genomic_prefixes = ('NC_', 'NT_', 'NW_', 'AC_', 'NZ_')
    print("Parsing catalog: {}".format(catalog_file))
    catalog_dic = {}
    with optionally_compressed_handle(catalog_file, 'rb') as fin:
        for line in fin:
            fields = line.split('\t')
            taxid = int(fields[0].strip())
            accession = str(fields[2].strip())
            if len(fields) == 6 and accession.startswith(genomic_prefixes):
                size = int(fields[5].strip())
                catalog_dic[accession] = (taxid, size)
            elif len(fields) == 7 and accession.startswith(genomic_prefixes):
                size = int(fields[6].strip())
                catalog_dic[accession] = (taxid, size)
    return catalog_dic

# SWITCH TO SQL
def my_awesome_func(output_dir, target_version):
    current_version = get_current_release_version(output_dir)
    # Always get the current catalog
    current_catalog = get_refseq_release_catalog(output_dir, current_version)
    current_dic = genomic_records_to_dic(current_catalog)

    if current_version != target_version:
        target_catalog = get_refseq_release_catalog(output_dir, target_version)
        target_dic = genomic_records_to_dic(target_catalog)
    else:
        target_dic=current_dic

    current_accession_set, target_accession_set = set(current_dic.keys()), set(target_dic.keys())
    common_accessions = target_accession_set.intersection(current_accession_set)
    unique_to_target = target_accession_set.difference(current_accession_set)
    unique_to_current = current_accession_set.difference(target_accession_set)
    print(map(len, [common_accessions, unique_to_current, unique_to_target]))


# get current release catalog
# is it the target?
# if not get also the target
# reconstitute the target catalog if necessary
# emit the sets unique to current, unique to previous and overlapping