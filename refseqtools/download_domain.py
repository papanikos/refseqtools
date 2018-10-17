import logging
import os
from .utils import NcbiFTPConnector, check_or_create_dir, select_genomic_fnas
from multiprocessing import Pool
from functools import partial
import hashlib

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def get_md5sum(filename):
    with open(filename, 'rb') as fin:
        file_as_bytes = fin.read()
    md5sum = hashlib.md5(file_as_bytes).hexdigest()
    return md5sum


def write_fastas_md5(input_dir):
    md5sums_file = os.path.join(input_dir, 'md5sums.txt')
    md5sums_dic = {}
    for filename in os.listdir(input_dir):
        if filename.endswith('fna.gz'):
            abs_fp = os.path.join(input_dir, filename)
            file_md5 = get_md5sum(abs_fp)
            md5sums_dic.update({abs_fp: file_md5})

    with open(md5sums_file, 'w') as fout:
        for f, md5 in md5sums_dic.items():
            fout.write('{}\t{}\n'.format(f, md5))

# TO DO: Revisit the whole function.
# It works when a stable connection is available.
# The actual process of downloading is not monitored,
# so that restarts when and if needed can be done.
def download_domain_file_list(base_dir, domain_name, source_files):
    """
    Given a list of files, download them from the NCBI ftp.
    Can also perform optional filtering after the raw_file has been written.

    :param base_dir: The directory where the files will be stored
    :param domain_name: The directory where the files will be downloaded from
    :param source_files: List of files to be downloaded.
    :return:
    """
    # Make a new connection for the process
    conn = NcbiFTPConnector()
    # Change cwd to the specified domain
    conn.go_to_dir(domain_name)

    # # Check or create the output directory
    # logger.info("Checking if {} exists".format(base_dir))
    target_dir = os.path.join(base_dir, domain_name)
    check_or_create_dir(target_dir)

    # This downloads the files in the list
    for source_file in source_files:
        target_file = os.path.join(target_dir, source_file)
        logger.info("Writing file {}".format(target_file))
        with open(target_file, 'wb') as fout:
            conn.ftp.retrbinary('RETR %s' % source_file, fout.write)

    # Close the ftp connection
    conn.ftp.quit()


def download(domain, output_dir, no_of_processes=4):

    # Make a connection to get the list of files to be downloaded
    c = NcbiFTPConnector()
    c.go_to_dir(domain)
    domain_files_list = c.ftp.nlst()
    c.ftp.quit()

    # # Check or create the output directory
    domain_dir = os.path.join(output_dir, domain)
    logger.info("Checking if {} exists".format(domain_dir))

    check_or_create_dir(domain_dir)

    # Select only the .genomic.fna.gz files
    genomic_fnas = select_genomic_fnas(domain_files_list)
    logger.info("{} files will be downloaded to {}".format(len(genomic_fnas), domain_dir))

    parallel_download = partial(download_domain_file_list, output_dir, domain)

    # PARALLEL DOWNLOADING

    # Create batches of 10 for parallel downloading
    genomic_fnas_batches = [genomic_fnas[i:i + 10]
                            for i in range(0, len(genomic_fnas), 10)]

    # If there are not more than 10 files, do not start multiple threads
    if len(genomic_fnas_batches) == 1:
        logging.info("Single process mode...")
        parallel_download(genomic_fnas_batches[0])
    else:
        # If the number of batches is smaller than the defined number of processes
        # then select the number of batches as the number or processes
        no_of_processes = min(len(genomic_fnas_batches), no_of_processes)

        logger.info("Starting {} parallel process for downloading..".format(no_of_processes))
        pool = Pool(processes=no_of_processes)
        pool.map(parallel_download, genomic_fnas_batches)
        pool.close()

    logger.info('Collecting md5sums for downloaded files')

    write_fastas_md5(domain_dir)
