#!/usr/bin/env nextflow

//Pipeline for flye, minimap, samtools, medaka, and augusutus
// first writing in one script, then breaking down into modules

//set number of cpus to reflect hardware via nextflow syntax
//update threads in each script accordingly

//input: *.fq.gz file



process flye {
    container 'quay.io/biocontainers/flye:2.9.5--py310h275bdba_2'
    input:
        path 'reads.fq.gz'
    output:
        path 'output/assembly.fasta', emit: fasta
        path 'output/flye.log'
        path 'output/assembly_info.txt'
    script:
    """
    flye --nano-hq reads.fq.gz  --out-dir ./output/ --genome-size 37m  --threads $params.cpus --asm-coverage 30
    """
}




process minimap {
    container 'quay.io/biocontainers/minimap2:2.17--hed695b0_3'

    input:
        path 'reads.fq.gz'
        path 'assembly.fasta'

    output:
        path 'minimap.bam'
    script:
    """
    minimap2 -ax map-ont -t $params.cpus assembly.fasta reads.fq.gz > minimap.bam
    """

}

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



process medaka {
    container 'quay.io/biocontainers/medaka:1.4.4--py38h130def0_0'
    input:
        path 'reads_to_draft.bam'
        path 'reads_to_draft.bam.bai'
        path 'assembly.fasta'
    output:
        path 'assembly_polished/consensus.fasta'

    script:
    """
    medaka_consensus -i reads_to_draft.bam -d assembly.fasta -o assembly_polished -t $params.cpus
    """

}


// more lines of Augustus code to get more files with extra content if we want to expand the pipeline
process augustus {
    container 'quay.io/biocontainers/augustus:3.5.0--pl5321heb9362c_5'
    input:
    path 'consensus.fasta'

    output:
    path 'gene_calling.gff3'

    script:
    """
    augustus --species=aspergillus_oryzae --strand=both --genemodel=partial --gff3=on consensus.fasta > gene_calling.gff3
    """

}


workflow{
    reads= file(params.input)
    flye(reads)
    minimap(reads,flye.out.fasta)
    samtools(minimap.out)
    medaka(samtools.out.bam,samtools.out.bai,flye.out.fasta)
    augustus(medaka.out)
}
