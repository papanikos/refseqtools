#!/usr/bin/env python

import sqlite3
from utils import *
import argparse

parser = argparse.ArgumentParser(description="Create an sqlite.db from a processed refseq catalog file"
                                 )
parser.add_argument('-c',
                    dest='catalog',
                    help="Path to a processed refseq release catalog, containing 3 columns: "
                         "accession, taxid, size",
                    required=True
                    )
parser.add_argument('-db',
                    dest='db_path',
                    help='Full path to output sqlite.db file',
                    required=True
                    )

if __name__ == '__main__':
    args = parser.parse_args()
    sqlite_file = args.db_path
    # Define the name of the table
    table_name = "acc_sizes"

    # Initialize a connection to the sqlite file
    conn = sqlite3.connect(sqlite_file)
    c = conn.cursor()

    # Create the table with an sql query
    c.execute("CREATE TABLE IF NOT EXISTS {tn} "
              "(acc text PRIMARY KEY,taxid integer,size integer);".format(tn=table_name))

    # Populate the rows of the table with the
    # refseq catalog data
    with optionally_compressed_handle(args.catalog, 'r') as fd:
        for l in fd:
            fields = l.split("\t")
            acc = fields[0].strip()
            taxid = int(fields[1].strip())
            size = int(fields[2].strip())
            c.execute("INSERT INTO {tn}(acc, taxid, size) VALUES(?,?,?);".format(tn=table_name), (acc, taxid, size,))

    conn.commit()
    conn.close()
