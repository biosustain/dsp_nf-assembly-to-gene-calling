# dsp_nf-assembly-to-gene-calling

This nextflow pipeline can take sequencing data from Oxford Nanopore Technologies (fastq or gzipped fastq) and conduct required steps in order to attain gene prediction as a gff3 file.

The pipeline relies on the following software that have been integrated in this nextflow pipeline:
* Flye - _de novo_ assembler for single-molecule sequencing reads - https://github.com/mikolmogorov/Flye
* minimap2 - Aligning reads to assembly - https://github.com/lh3/minimap2
* samtools - Sorting and indexing BAM - https://github.com/samtools/samtools
* medaka - Polishing (consensus sequence) - https://github.com/nanoporetech/medaka
* augustus - Gene prediction - https://github.com/Gaius-Augustus/Augustus

 currently only takes the following input parameters:

* input: *.fastq.gz file
* threads: integer - how many threads to parallelize processes that support it
* outdir -  to be used as publishDir
* species: species to be used by Augustus (lower case, underscore instead of space)
* genome_size - for Flye
* assembly_coverage - for Flye


example to initiate run:

```
nextflow run https://github.com/biosustain/dsp_nf-assembly-to-gene-calling --input my_data.fq.gz --species aspergillus_oryzae --genome_size 37m --assembly_coverage 30 --threads 16 -resume
```
