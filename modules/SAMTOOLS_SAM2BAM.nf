process SAMTOOLS_SAM2BAM {
    label 'midmem'
    container 'ebird013/samtools:1.17'

    input:
        tuple val(sample), file(alignment)
    output:
        tuple val(sample), path("${sample}.bam"), emit: bam

    script:

    """
    samtools view -bS ${alignment} > ${sample}.bam
    """
}