#!/usr/bin/env python

import argparse
import logging
from utils import *
from multiprocessing import Pool
from functools import partial

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser(description=
                                 "Download genomic.fna.gz files "
                                 "from the ncbi ftp RefSeq release for a given domain")
parser.add_argument('-d',
                    dest='domain',
                    help="One of 'viral', 'bacteria', 'invertebrate', 'archaea' "
                         "'fungi', invertebrate', 'protozoa', 'vertebrate_mammalian', "
                         "'vertebrate_other'",
                    required=True
                    )
parser.add_argument('-o',
                    dest='output_dir',
                    help='Full path to output directory',
                    required=True
                    )
parser.add_argument('-p',
                    dest='no_of_processes',
                    type=int,
                    help='Number of processes to start.',
                    default=4,
                    required=False)


def download_domain(base_dir, domain_name, source_files):
    """
    Given a list of files, download them from the NCBI ftp.
    Can also perform optional filtering after the raw_file has been written.

    :param base_dir: The directory where the files
    :param domain_name:
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


def main():
    # Parse the arguments
    args = parser.parse_args()

    # Make a connection to get the list of files to be downloaded
    c = NcbiFTPConnector()
    c.go_to_dir(args.domain)
    domain_files_list = c.ftp.nlst()
    c.ftp.quit()

    # # Check or create the output directory
    domain_dir = os.path.join(args.output_dir, args.domain)
    logger.info("Checking if {} exists".format(domain_dir))

    check_or_create_dir(domain_dir)

    # Select only the .genomic.fna.gz files
    genomic_fnas = select_genomic_fnas(domain_files_list)
    logger.info("{} files will be downloaded to {}".format(len(genomic_fnas), domain_dir))

    parallel_download = partial(download_domain, args.output_dir, args.domain)

    ## PARALLEL DOWNLOADING

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
        no_of_processes = min(len(genomic_fnas_batches), args.no_of_processes)

        logger.info("Starting {} parallel process for downloading..".format(no_of_processes))
        pool = Pool(processes=no_of_processes)
        pool.map(parallel_download, genomic_fnas_batches)
        pool.close()


if __name__ == '__main__':
    main()
