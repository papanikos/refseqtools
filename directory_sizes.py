#!/usr/bin/env python

from utils import NcbiFTPConnector
from collections import namedtuple

# Domains of interest
DOMAINS = ['archaea', 'bacteria', 'complete', 'fungi', 'invertebrate',
           'mitochondrion', 'plant', 'plasmid', 'plastid', 'protozoa',
           'vertebrate_mammalian', 'vertebrate_other', 'viral']


def get_files_sizes(domain):
    """
    Get the file sizes for a specified domain from /refseq/release/domain

    :param domain: Domain for which the stats will be retrieved
    :return: A named tuple (FileInfo(domain, no_of_files, total_size
    """
    # Initiate a connection to the NCBI ftp site
    conn = NcbiFTPConnector()
    conn.go_to_dir(domain)
    # Store the file information in a list
    dir_info = []
    conn.ftp.retrlines('LIST', dir_info.append)
    # Close the connection
    conn.ftp.quit()

    # Initialize counters
    no_of_files = 0
    total_size = 0

    # Loop over the files list and get the information
    for l in dir_info:
        fields = l.split()
        filename = fields[-1]
        filesize = int(fields[4])
        if filename.endswith('genomic.fna.gz'):
            no_of_files += 1
            total_size += filesize

    # Calculate the total size in GB
    total_size_gb = round(total_size/(1024**3), 3)

    # Instantiate the named tuple
    FileInfo = namedtuple('FileInfo', ['domain', 'no_of_files', 'total_size'])
    domain_info = FileInfo(domain=domain,
                           no_of_files=no_of_files,
                           total_size=total_size_gb)
    return domain_info


if __name__ == '__main__':
    print("#{:<20}\t{: >10}\t{: >10}".format("Domain", "No. of files", "Size (compressed, GB)"))
    for domain in DOMAINS:
        a = get_files_sizes(domain)
        print('{:<20}\t{: >10}\t{: >10}'.format(a.domain, a.no_of_files, a.total_size))
