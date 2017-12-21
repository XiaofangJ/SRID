# SRID ("Split Read Insertion Detection”)

## Overview
SRID is a pipeline wrapped in python that uses samtools, bedtools and BEDOPS to 
identify putative inserted regions (candidate mobile genetic elements) in 
reference genomes and extracts reads supporting such insertions. 
SRID is the first step to Identify MGEs in reference genomes. 

If you have any questions about SRID, feel free to contact me (xiaofang AT mit DOT edu).

The algorithm implemented in this pipeline is similar to 
[Daisy](https://github.com/ktrappe/daisy.git). To my knowledge, Daisy is the 
first tool to perform  horizontal gene transfer detection by mapping sequencing reads.

## Software requirement 
+ [samtools](http://samtools.sourceforge.net/) (>=1.4)
+ [bedtools](http://bedtools.readthedocs.io/en/latest/) (>=v2.26.0)
+ [BEDOPS](https://bedops.readthedocs.io/en/latest/) (>=2.4.20)


## Quick Start
```
usage: SRID.py [-h] -b  -o  [-p] -r  -m  -s  [-n]

Identify putative insertions in reference genomes from bam alignments

optional arguments:
  -h, --help       show this help message and exit
  -b , --bam       input bam file
  -o , --out       output file
  -p , --threads   number of threads used
  -r , --readlen   average read length
  -m , --mean      mean library insertion size
  -s , --sd        standard deviation of library insertion size
  -n , --number    number of split reads required to validate split sites
```

### Input
+ input bam file: bam alignments generated from 
`bwa mem ref.fa 1.fq 2.fq |samtools view -F 4 -Sb >input.bam` 
     + ref.fa ==> reference genomes where MGEs are expected to inserted into 
     + 1.fq, 2.fq ==> paired files with metagenomic/genomic sequencing read pairs;
     the name of the reads should not be identical in 1.fq and 2.fq
+ average read length: the mean value of read length in 1.fq and 2.fq 
+ library insertion size: the base pairs of the sequence between adapters 
    + You can either obtain insertion size information from your sequencing
    service provider or estimate it with tools such as:
       * CollectInsertSizeMetrics from [Picard](http://broadinstitute.github.io/picard/)
       * BAM QC from [QualiMap](http://qualimap.bioinfo.cipf.es/)

### Output
Output is a table of five columns.

 Column name           | Explanation                                                                                                                 |
-----------------------|-----------------------------------------------------------------------------------------------------------------------------|
Scaffolds              | The scaffolds where the insertions are detected
Split Site 1           | The start coordinates of insertions; The start of putative MGEs
Split Site 2           | The end coordinates of insertions; The end of putative MGEs
Supporting Split Reads | The list of reads separated by commas that split and are aligned across the insertion sites
Supporting Read Pairs  | The list of read pairs (only show the name of the first in pair) separated by commas that flank and support the insertions.

Example: 
```
python SRID.py -b ./data/example/test.bam -p 12 -r 100 -m 200 -s 71 -n 4 -o out.tab
```

The output should be the same as `./data/example/out.tab` which has contains two entry :

NZ_JH724283.1:750858-802013 is an ICE        [check it](https://immedb.mit.edu/NZJH7242831750858-802013) 
NZ_JH724283.1:844019-845389 is a transposon   [check it](https://immedb.mit.edu/NZJH7242831844019-845389)  


## Notes
1. Daisy integrates coverage and read pair information for HGT candidate evaluation. 
It is a tool intended to detect HGT with high confidence. You should use 
Daisy if your sequencing data is from single organism with suspected HGT 
events and you have reference genomes for acceptor and donor genomes.  
2. SRID is designed to detect putative insertion in reference genomes. SRID 
doesn't use coverage information to verify putative insertion region, as the 
coverage will be highly variable in metagenomic context. 
3. SRID is intended to identify as many as putative MGEs as possible and especially
in low-abundance species in metagenomic sequencing data. Therefore, there are no 
stringent requirements for filtering based on read or mapping quality.  
You should always check the reliability of the putative region with the following methods: 
    + Extract those reads supporting such putative regions from the output file 
    to check their quality by means such as visualizing the alignment.
         + [Extract reads from bam file by read name](http://timoast.github.io/2015/10/12/ExtractReads/)
    + Search for MGE gene signatures to validate and classify the putative MGEs

## Citation 
Xiaofang Jiang, Andrew Brantley Hall, Ramnik J Xavier, Eric J Alm (2016)
Comprehensive analysis of mobile genetic elements in the gut microbiome reveals a phylum-level niche-adaptive gene pool
bioRxiv 214213; doi: https://doi.org/10.1101/214213i
