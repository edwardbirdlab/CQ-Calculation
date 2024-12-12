process SAMTOOLS_COUNTS {
    label 'verylow'
    container 'ebird013/samtools:1.17'

    input:
        tuple val(sample), file(alignment)
    output:
        tuple val(sample), path("${sample}_counts.txt"), emit: stats

    script:

    """
    samtools idxstats ${alignment} > ${sample}_counts.txt
    """
}