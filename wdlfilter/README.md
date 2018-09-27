# WDL workflows for filtering fasta files

`FilterDomain.wdl` is a generalized workflow (depending on the FilterDomainTasks.wdl) that can be used to scatter the
filtering jobs and the dustmasking of all the previously downloaded data files.

An example for filtering fasta files downloaded for the domain of invertebrates
```
$ java -jar cromwell.34.jar run FilterDomain.wdl -i invertebrate.json
```

## Filtering bacteria

Most RefSeq genomic sequences for bacteria are NZ accessions, resulting from WGS experiments.
Based on recommendations from the RefSeq maintainers, a list of representative assemblies (and the accessions included
in these) can be retrived with the following query on the NCBI database:

https://www.ncbi.nlm.nih.gov/assembly/?term=prokaryotes%5BOrganism%5D+AND+%22latest+refseq%22%5Bfilter%5D+AND+(%22representative+genome%22%5Bfilter%5D+OR+%22reference+genome%22%5Bfilter%5D)

The assemblies reports that include the accessions can be retrieved. A list of the accessions
can be created and used as input to the filter_fasta.py script, in order to retain these sequences.
In addition, accessions corresponding to other chromosomes and plasmids (prefixed `NZ_CM` and `NZ_CP`) were retained.


## Redundant accessions

Many NC (chromosome) accessions are made up from several NW and NT (scaffolds) accessions. These sequences are all
included redundantly in the downloaded fasta files. A list of these accessions can be retrieved with this query on the
NCBI databases

https://www.ncbi.nlm.nih.gov/nuccore/?term=(NT_000001%3ANT_999999%5Bpacc%5D+OR+NW_000001%3ANW_999999%5Bpacc%5D)+AND+nuccore_comp_nuccore%5Bfilter%5D

These NW and NT accessions are distributed across different domains (invertebrate, plants, vertebrate_other and 
vertebrate_mammalian)as shown [here](https://barmsijs.lumc.nl/Nikos/refseq_pics/nuccore_result.png).

The aforementioned list was passed as input to the filter_fasta.py script for each of these domains.

## Filtering vertebrate-mammalian

We decided to keep RefSeq records belonging only to certain species, namely:
 - _H. sapiens_ (human, NCBI taxid: 9606)
 - _M. musculus_ (mouse, NCBI taxid: 10090)
 - _R. norvegicus_ (rat, NCBI taxid: 10116)

that can be relevant in the clinical context.


