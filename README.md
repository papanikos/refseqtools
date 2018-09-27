# download_refseq_release

Scripts for handling NCBI taxonomy and RefSeq Release information.

# Environment

For [handling the NCBI taxonomy](http://etetoolkit.org/docs/latest/tutorial/tutorial_ncbitaxonomy.html)
 (and [tree processing](http://etetoolkit.org/docs/latest/tutorial/index.html) thereof)
  the [ETE3 toolkit](http://etetoolkit.org/docs/latest/index.html)
  is used

```
conda env create -f=environment.yml
```

Warning: the ete3 toolkit and its dependencies take up almost 500MB.

All subsequent commands assume you are in the activated environment `source activate ete3`

> Note on using ete3's NCBITaxa()
>
>When creating a tree object, ete3's NCBITaxa() is called. When this is done for the first time, taxonomy related 
information is downloaded from scratch by default in the home directory. 
>If you already have an instance of ete3's created sqlite DB
>in another directory you can use that with the `-db` option.

## Refseq Release catalog

Most scripts rely on a json file that contains sequence information for each taxonomy available from
[here](ftp://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/).
The json file looks like this:
```
{taxid : { accession: size, ... }, ... }
```

This was created with:
```
# Download the file
$ wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/RefSeq-release89.catalog.gz

# Filter only DNA refseq sequences
$ zgrep -e "AC_" -e "NC_" -e "NG_" -e "NT_" -e "NW_" -e "NZ_" | \
awk -F "\t" {print $3,"\t",$1,"\t",$6} | \
gzip -c > Refseq.89.DNA.gz

# Create the file with the python script
$ python convert_acc2taxid_to_json.py -i Refseq.89.DNA.gz -o Refseq.89.json

```

> Improvement 
>
>`convert_acc2taxid_to_json.py` can be modified to skip the intermediate
>preprocessing step and retain accessions and taxids according to the
>user's needs 


# Getting domain sizes - directory_sizes.py

- This will get the directory sizes of the current refseq release.
- It can serve as an estimate of the total download size for each domain.
- `directory_sizes.py` prints to stdout

```
## RefSeq v89
$ python directory_sizes.py

#Domain                 No. of files    Size (compressed, GB)
archaea                          6           0.593
bacteria                      1377         135.105
complete                      2835         341.798
fungi                           59           2.198
invertebrate                   149           17.72
mitochondrion                    2           0.062
plant                          178          28.068
plasmid                          4           0.367
plastid                          2           0.125
protozoa                        19           0.986
vertebrate_mammalian           715          104.06
vertebrate_other               313          52.983
viral                            2           0.077

```

```
## RefSeq v90
$ python directory_sizes.py

#Domain                 No. of files    Size (compressed, GB)
archaea                          7           0.622
bacteria                      1484         150.773
complete                      2978         363.974
fungi                           60           2.283
invertebrate                   154          17.735
mitochondrion                    2           0.064
plant                          190          29.739
plasmid                          5           0.386
plastid                          2           0.129
protozoa                        21             1.0
vertebrate_mammalian           733         106.121
vertebrate_other               327          55.614
viral                            3           0.078
```

# Download sequences for a domain - download_refseq_domain.py

```
$ python download_refseq_domain.py -d viral -o path/to/output_dir
```

will download all viral.*.genomic.fna.gz files to `output_dir/viral`. 

The `-p` option can be used to start multiple download processes in parallel.

# Annotate NCBI taxonomy tree with sequence information - annotate_tree.py

- This requires only the information included in the release catalog - i.e. any downloaded fasta files are not 
required.

To get an annotated tree with sequence information for Hominidae (taxid= 9604) pruned at the species level.
```
$ python annotate_tree.py -t 9604 -j path/to/Refseq.json -p species
```

To write the image to a `.png` file (supported image output options are `.png`, `.svg`, `.pdf`)
```
$ python annotate_tree.py -t 9604 -j path/to/Refseq.json -p species -r path/to/image.png
```

Writing a newick tree for use with other tree handling software (e.g. [iTOL](https://itol.embl.de/)) is also possible,
along with the `render` option. The output newick tree will have the same basename as the rendered image, suffixed
with `.nw`.

```
$ python annotate_tree.py -t 9604 -j path/to/Refseq.json -r path/to/9604.png -write-to-file

# The newick formatted tree will be written to path/to/9604.nw
```

# Filter a RefSeq fasta file - filter_fasta.py

- Enables filtering of an input fasta file based on:
  - Accession type prefix (e.g. NC, NZ) of the fasta records
  - Taxonomy - keep/remove fasta records belonging to user defined taxonomies
  - Predefined lists of accessions
  
For example, if you would like to filter a fasta file containing viral genomic sequence records, excluding all
accessions prefixed with `AC_` and including records belonging to Adenoviridae (taxid: 10508) or Orthomyxoviridae
(taxid: 11308):
  
```
$ python filter_fasta -i path/to/viral.1.1.genomic.fna.gz \
--exclude-accessions AC \
--include-taxa 10508,11308
-j path/to/Refseq.json
```

# Switching to a sqlite db

The Refseq.json that contains the sequence information per accession might not be the best
option to store this information. As an alternative, the refseq release catalog can be
stored in a sqlite database, that is more flexible and performs better, instead of having
to load all the information each time.

This option for now is explored with the `annotate_tree_from_fasta.py`, where given a fasta
file containing Refseq genomic records a tree representing the contents of the fasta file 
can be drawn.

To achieve this, first create the database

```
$ python create_refseq_db.py -c path/to/refseq_catalog -db path/to/sqlite.db
```

>Note
>
>The input for this script is assumed to be a processed Refseq release catalog, containing
>3 columns: accession, taxid, size.

After the database has been successfully created, an annotated tree with the sequence
content for each node can be drawn, similar to the `annotate_tree.py` module.

```
$ python annotate_tree_from_fasta.py -i path/to/input.fasta \
-acc_db path/to/refseq_release.db \
-r path/to/image.png
```

>Improvement
>
>Ideally, all other scripts should be refactored to use this format.

# Plotting redundant accessions

As mentioned [here](wdlfilter/README.md#Redundant-accessions), several chromosome accessions
are in fact joins of scaffold accessions. In order to explore which domains these are
included in, the `plot_nuccore_summary.py` was used. 

The input file for this is a summary.txt file, downloaded after running the relevant query.

