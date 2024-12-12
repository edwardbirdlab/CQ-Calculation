process SAMTOOLS_COUNTS {
    label 'verylow'
    container 'ebird013/samtools:1.17'

    input:
        tuple val(sample), file(alignment)
    output:
        tuple val(sample), path("${sample}_align.txt"), emit: stats

    script:

    """
    samtools view -Sb ${alignment} > ${sample}.bam
    samtools index ${sample}.bam
    samtools idxstats ${sample}.bam> ${sanple}_counts.txt
    """
}