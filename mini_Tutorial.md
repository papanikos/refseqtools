# Introduction

The idea is to use standardized NCBI information (the RefSeq Release source of
sequence information, NCBI taxonomy) to ensure reproducibility of the index 
building process, for use with metagenomics classifiers.

# Building a centrifuge index

`centrifuge-build` requires as input:
- A set of input sequences
- A taxonomy tree
- A mapping of sequences to their corresponding organism in the tree
- A mapping of the taxonomy ids to meaningful names

## Requirements
The information covered here referes to the old_master branch.
To get the necessary scripts and workflows you should run
`git clone --single-branch -b old_master https://github.com/papanikos/refseqtools.git`.
This will also fetch the `requirements.yml` required to set up a conda environment.


# An example

We will try to build a centrifuge index for all viruses,
 based on the current refseq release, filtered based
on the following criteria:
- Exclude all sequences belonging to the family of Adenoviridae (taxid: 10508)
- Exclude all `AC_` accessions


## Step 1a. Getting the refseq-release catalog

A catalog that summarizes all the contents of the release is available from NCBI's ftp site 
ftp://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/

Retrieving the catalog
```
$ wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/RefSeq-release90.catalog.gz
```

This contains all the information for the sequences included in the current release.
This includes reference sequences for genomic records, but also proteins, mRNAs.
In this context we are not really interested in all of the records, but we would
like to focus on only the genomic ones.

Refseq accessions are logically named following certain conventions. These are
summarized in [this table](https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/?report=objectonly).

We can see that some of the accessions can actually be excluded for further
analyses. By reducing the size of the catalog, not only do we retain the relevant
information for our goal, but we also reduce the size of I/O.

The following chain of commands will grab only the genomic accessions, retain only
the accession (column 3), its taxid (column 1) and its size (column 6) and rearrange
the output in a new gzipped table.

```
$ zgrep -e "AC_" -e "NC_" -e "NW_" -e "NT_" -e "NZ_" RefSeq-release90.catalog.gz | \
awk -F "\t" '{print $3"\t"$1"\t"$6}' | gzip -c > Refseq90.DNA.gz
```

The new table now contains the filtered input in rearranged columns
```
$ zcat Refseq90.DNA.gz | head -n 20
NC_001911.1	9	7768
NC_004843.1	9	2308
NC_013549.1	9	6521
NC_025017.1	9	3037
NZ_LT667500.1	9	447673
NZ_LT667501.1	9	6126
NZ_LT667502.1	9	2904
NZ_MCBT00000000.1	23	4642463
NZ_MCBT01000001.1	23	325457
NZ_MCBT01000002.1	23	6636
```

From this table we can already extract the mapping file required from centrifuge
that contains a mapping of each sequence to its taxid.

```
zcat Refseq90.DNA.gz | cut -f1,2 > seqid2taxid.map 
```

Let's also put this into a more convenient data structure (a json file)
that can be readily loaded into scripts that use this information to calculate
sequence contents per taxonomy. This is done with the `convert_acc2taxid_to_json.py`
utility, where a taxid is used as a primary key and all accessions and sequence
sizes are stored under it.

```
$ python convert_acc2taxid_to_json.py -i Refseq90.DNA.gz -o Refseq90.DNA.json 
```

## Step 1b. Getting the taxonomy tree

Now we have the information about the taxonomy of each genomic sequence and its size
conveniently stored in one place. But we still need a taxonomy tree, that
describes the relationships between the taxids.

For this we are using the ete3 package, which provides:
- an interface to the NCBI taxonomy
- tools for manipulating phylogenetic trees for visualization

First we need to instantiate a SQL DB that ete3 uses for the NCBI taxonomy
manipulation.

  1. Retrieve the taxonomy information from the NCBI taxonomy DB, located
here ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz .

This is a tarball that contains all necessary information.

>The NCBI taxonomy is regularly updated. Downloading this at different time points
> will affect the end-result of the database creation and all analyses thereafter.

```
$ wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
```

We can also extract the required `names.dmp` and `nodes.dmp` for later use.

```
tar -zxvf taxdump.tar.gz names.dmp nodes.dmp
```

  2. Create the database
 
ete3 requires a sqlite db file to already exist. We need to first create this
(empty) database and then point ete3 to the location of the taxdump file (
downloaded from above) we want
it to use to populate the database. For demonstration purposes, this is done here
with an interactive python session
```
(ete3)$ python
Python 3.6.6 | packaged by conda-forge | (default, Jul 26 2018, 09:53:17) 
[GCC 4.8.2 20140120 (Red Hat 4.8.2-15)] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> import sqlite3
>>> from ete3 import NCBITaxa
>>> connection = sqlite3.connect("/path/to/taxa.sqlite") # Change the paths
>>> connection.close()
>>> NCBITaxa(dbfile="/path/to/taxa.sqlite", taxdump_file="/path/to/taxdump.tar.gz")
Loading node names...
1934824 names loaded.
205002 synonyms loaded.
Loading nodes...
1934824 nodes loaded.
Linking nodes...
Tree is loaded.
Updating database: taxa.sqlite ...
 1934000 generating entries... 
Uploading to taxa.sqlite

Inserting synonyms:      205000 
Inserting taxid merges:  50000 
Inserting taxids:       1930000 
```

Now we have all required inputs for manipulating the taxonomy and the
sequences accessions in order to build our index.

## Step 2. Summarizing information and plotting

Before we start downloading the actual sequences it might be helpful to get
an idea of what the release contents for our chosen taxonomies actually are.

Since we are interested in creating an index based on the viral sequences only,
let's have a closer look.

The `annotate_tree.py` script provides this functionality. For now it relies on
taxids, so we need to know what the actual taxid for the whole group of viruses
is. A quick search on the NCBI taxonomy site shows that this is 10239.

If we now feed this into the script we can create a tree visualization that 
will summarize all the information for us

> Remember, we created the necessary input in the previous steps.
> These would be the Refseq90.DNA.json and the NCBI taxonomy database.
> The full paths to their locations are required.

```
$ python annotate_tree.py -t 10239 \
-db /path/to/ete3_db/ \
-j /path/to/Refseq90.DNA.json \
-r viruses.svg
```

The created `viruses.svg` can be opened within a browser window. It summarizes
all the sequences contained under the Viruses taxonomy, up to the species level.
This is a very detailed view and somewhat hard to completely browse through. We can choose
the level at which the same information is visualized with the `-prune` option
of the script.

``` 
$ python annotate_tree.py -t 10239 \
-db /path/to/ete3_db/ \
-j /path/to/Refseq90.DNA.json \
-r viruses.family.svg \
-prune family 
```

## Step 3. Downloading sequences

We had a view of the total information that is contained within the kingdom of
viruses. We can now retrieve the fasta files these are contained in.

This is done with the `download_domain.py` utility.

```
$ python download_refseq_domain.py -d viral -o /path/to/output_dir/
```

Once the download is complete, there should be 3 fasta files (for release 90) including all the
reference sequences for the domain of viruses.

## Step 4. Filtering the files according to our needs

We set out to make a reference set that should include all viral sequences,
except for Adenoviridae and AC_ accessions.
The taxid for Adenoviridae is 10508. 

- For one file

The `filter_fasta.py` enables us to filter the sequences from a selected fasta
file, based on the desired criteria.

>Important
>
>In its current implementation, the script modifies the selected file.
>If you want to retain the original files intact, you should copy them
>some place safe and work with copies of the files.

If we would like to filter one of the downloaded files according to our criteria
```
$ python filter_fasta.py \
-i /path/to/downloaded/viral.2.1.fna.gz \
-o /path/to/output/filtered_fasta.fna.gz \
--exclude-taxa 10508 \
--exclude-accessions AC \
-j /path/to/Refseq90.DNA.json   
```

- For multiple files

When dealing with thousands of files over different domains, running the script
for each file is facilitated with wdl defined workflows, located in the 
[wdlfilter directory](./wdlfilter).

The workflow FilterDomain.wdl assumes that a domain directory is given, where all
the desired fna files, downloaded from refseq release are in. Then applies the
filtering criteria to retain/remove sequences from the fasta files.

It further expands on the filtering to add dustmasking on the resulting files and
concatenates all files. The resulting file should be ready to be input to centrifuge
for building the customized viral index.

>In its current status, the workflow runs under cromwell v34.
>The dustmasker executable also needs to be installed independently.

For running the workflow an `inputs.json` file is required which 
configures the pipeline. To produce such a file, the womtool-34 jar is required,
available [here](https://github.com/broadinstitute/cromwell/releases/tag/34).
The cromwell-34.jar is also required to run the complete workflow, once the inputs
are defined.

Once the jars are downloaded, head over to the downloaded wdlfilter directory and run
```
cd wdlfilter/
java -jar /path/to/womtool-34.jar inputs FilterDomain.wdl > viral.inputs.json 
```

The resulting file should look like
``` 
{
  "FilterDomainFastas.filterScriptPath": "String",
  "FilterDomainFastas.Cat.zip": "Boolean (optional, default = false)",
  "FilterDomainFastas.includeAccFile": "File? (optional)",
  "FilterDomainFastas.domainInputDir": "String",
  "FilterDomainFastas.excludeAccFile": "File? (optional)",
  "FilterDomainFastas.taxDbPath": "String",
  "FilterDomainFastas.includeAccPrefixes": "Array[String]? (optional)",
  "FilterDomainFastas.excludeAccPrefixes": "Array[String]? (optional)",
  "FilterDomainFastas.FilterFasta.threads": "Int (optional, default = 1)",
  "FilterDomainFastas.excludeTaxa": "Array[String]+? (optional)",
  "FilterDomainFastas.outputDir": "String",
  "FilterDomainFastas.Dustmasker.isGzipped": "Boolean? (optional, default = true)",
  "FilterDomainFastas.FilterFasta.preCommand": "String? (optional)",
  "FilterDomainFastas.Cat.unzip": "Boolean (optional, default = false)",
  "FilterDomainFastas.includeTaxa": "Array[String]+? (optional)",
  "FilterDomainFastas.dustmaskerExe": "String",
  "FilterDomainFastas.refseqJson": "String"
  "FilterDomainFastas.dustmaskOnly: "Boolean (optional, default=false)"
}

```
These are all the configurable options. Options having an (optional) tag,
can be omitted but the rest should be filled in according to your setup. Note
that the options hint to the type of input that is expected.

A minimal working config file, for our case, would look like this (replace
the paths according to your directories)
```
{
  "FilterDomainFastas.filterScriptPath": "/path/to/filter_fasta.py",
  "FilterDomainFastas.domainInputDir": "/path/to/downloaded/domain/",
  "FilterDomainFastas.taxDbPath": "/path/to/ete3_db/taxa.sqlite",
  "FilterDomainFastas.excludeAccPrefixes": ["AC"],
  "FilterDomainFastas.excludeTaxa": ["10508"],
  "FilterDomainFastas.outputDir": "/path/to/filtered_fastas/",
  "FilterDomainFastas.FilterFasta.preCommand": "source activate ete3",
  "FilterDomainFastas.dustmaskerExe": "/path/to/ncbi-tools/bin/dustmasker",
  "FilterDomainFastas.refseqJson": "/path/to/Refseq90.json"
}
```

We can now run the workflow, to filter the fasta files, dustmask them and merge
the sequence files together.

```
java -jar /path/to/cromwell-34.jar run FilterDomain.wdl -i inputs.json 
```

If the workflow run succeeds, there should be a `dustmasked.filtered.fna.gz`
file within the input directory.
In case you would like to skip the filtering step and only dustmask the files
before running the index-building step you can set the `dustmaskOnly` option to `true`.

Since centrifuge-build expects an unzipped file
``` 
cd ../
zcat /path/to/filtered_fastas/dustmasked.filtered.fna.gz > /path/to/filtered_fastas/viral_sequences.fna
```

## Step 5. Building the index with centrifuge

We should now have all the information required by centrifuge to build the
index.

1. The taxonomy tree [`nodes.dmp` (step 1b)]
2. Descriptive names for the organisms [`names.dmp` (step 1b)]
3. A mapping between sequences-taxid [`seqid2taxid.map` (step 1a)]
4. A set of filtered sequences [`viral_sequences.fna` (step 4)]

Running `centrifuge-build` should now be straightforward:
``` 
centrifuge-build --conversion-table /path/to/seqid2taxid.map \
--taxonomy-tree /path/to/nodes.dmp \
--name-table /path/to/names.dmp \
/path/to/viral_sequences.fna viral
```

This should produce the `viral.*.cf` files that constitute the index. Now
we are ready to use this as an index for classifying input reads with centrifuge.


## Notes on NCBI taxonomy and taxids

- NCBI assigns unique taxonomic identifiers to each organism.
(eg. modern humans have taxid 9606)
- Taxids are not very descriptive to read so a human readable form
is the actual species name _Homo sapiens_
- All life is artificially divided in taxonomic ranks, that
progressively group more similar organisms together. Starting
from the root of the tree (taxid: 1) there are several kingdoms,
phyla, orders, families, genera and species (and various subgroups). Each one is unique
and a different taxid can be assigned to it. So the full taxonomic
path, using the NCBI taxonomy, leading to _Homo sapiens_ (taxid: 9606) would be

```
1, 131567, 2759, 33154, 33208, 6072, 33213, 33511, 7711, 89593, 7742,
7776, 117570, 117571, 8287, 1338369, 32523, 32524, 40674, 32525, 9347,
1437010, 314146, 9443, 376913, 314293, 9526, 314295, 9604, 207598, 9605,
9606
``` 
