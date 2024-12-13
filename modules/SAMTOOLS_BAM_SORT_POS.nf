process SAMTOOLS_BAM_SORT_POS {
    label 'samtoolssort'
    container 'ebird013/samtools:1.17'

    input:
        tuple val(sample), file(alignment)
    output:
        tuple val(sample), path("${sample}_sorted.bam"), emit: sort

    script:

    """
    samtools sort -m ${task.memory.toGiga()}G -@ ${task.cpus} ${alignment} -o ${sample}_sorted.bam
    """
}