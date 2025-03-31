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
* cpus: integer - how many threads to paralellise processes that support it
