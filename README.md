# dsp_nf-assembly-to-gene-calling

Pipeline assembly to gene calling with:
* Flye
* minimap2
* samtools
* medaka
* augustus

testing set up with _aspergillus oryzae_,
 currently only takes the following input parameters:

* input: *.fastq.gz file
* threads: integer - how many threads to parallelize processes that support it
* outdir -  to be used as publishDir
* species: species to be used by Augustus (lower case, underscore instead of space)
* genome_size - for Flye
* assembly_coverage - for Flye
