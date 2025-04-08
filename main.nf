#!/usr/bin/env nextflow

//set number of cpus to reflect hardware via nextflow syntax
//update threads in each script accordingly

/* 
 * default parameters
 */
params.threads = 16 //could be set to task.cpus, but seems to default at 1
params.input = 'reads.fq.gz'
params.outdir = 'output'
params.species = 'aspergillus_oryzae'
params.genome_size = '37m' // 37 megabases for Aspergillus oryzae
params.assembly_coverage = 30


/*
 * de novo assembler for single-molecule sequencing reads
 * source: https://github.com/mikolmogorov/Flye
 */
process flye {
    publishDir "${params.outdir}/flye", mode: 'copy', overwrite: true
    container 'quay.io/biocontainers/flye:2.9.5--py310h275bdba_2'
    input:
        path 'reads.fq.gz'
    output:
        path 'output/assembly.fasta', emit: fasta
        path 'output/flye.log'
        path 'output/assembly_info.txt'
        path 'output/assembly_graph.gfa'
    script:
        """
        flye --nano-hq reads.fq.gz  --out-dir ./output/ --genome-size $params.genome_size  --threads $params.threads --asm-coverage $params.assembly_coverage
        """
}

/*
 * Align Reads to Assembly with Minimap2
 * source: https://github.com/lh3/minimap2
 */
process minimap {
    container 'quay.io/biocontainers/minimap2:2.17--hed695b0_3'
    input:
        path 'reads.fq.gz'
        path 'assembly.fasta'

    output:
        path 'minimap.bam' 
    script:
        """
        minimap2 -ax map-ont -t $params.threads assembly.fasta reads.fq.gz > minimap.bam
        """
}

/*
 * Sort and Index BAM with Samtools 
 * source: https://github.com/samtools/samtools
 */
process samtools {
    container 'quay.io/biocontainers/samtools:1.15.1--h1170115_0'
    input:
        path 'minimap.bam'
    output:
        path 'reads_to_draft.bam', emit: bam
        path 'reads_to_draft.bam.bai', emit: bai
    script:
        """
        samtools sort -o reads_to_draft.bam minimap.bam
        samtools index reads_to_draft.bam
        """
}

/*
 * Consensus Sequences with Medaka 
 * source: https://github.com/nanoporetech/medaka
 */
process medaka {
    publishDir "${params.outdir}/medaka", mode: 'copy', overwrite: true
    container 'quay.io/biocontainers/medaka:2.0.1--py311hfd2b166_0'
    input:
        path 'reads_to_draft.bam'
        path 'reads_to_draft.bam.bai'
        path 'assembly.fasta'
    output:
        path 'assembly_polished/consensus.fasta', emit: fasta
    script:
        """
        medaka_consensus -i reads_to_draft.bam -d assembly.fasta -o assembly_polished -t $params.threads
        """
}

/*
 * Gene Prediction with Augustus
 * source: https://github.com/Gaius-Augustus/Augustus
 */
process augustus {
    publishDir "${params.outdir}/augustus", mode: 'copy', overwrite: true
    container 'quay.io/biocontainers/augustus:3.5.0--pl5321heb9362c_5'
    input:
        path 'consensus.fasta'
    output:
        path 'gene_calling.gff3'
    script:
        """
        augustus --species=${params.species} --strand=both --genemodel=partial --gff3=on consensus.fasta > gene_calling.gff3
        """
}
// more lines of Augustus code to get more files with extra content if we want to expand the pipeline

workflow{
    reads= file(params.input)
    flye(reads)
    minimap(reads,flye.out.fasta)
    samtools(minimap.out)
    medaka(samtools.out.bam,samtools.out.bai,flye.out.fasta)
    augustus(medaka.out.fasta)
}

