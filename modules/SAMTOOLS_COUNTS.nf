process SAMTOOLS_COUNTS {
    label 'samtoolssort'
    container 'ebird013/samtools:1.17'

    input:
        tuple val(sample), file(alignment)
    output:
        tuple val(sample), path("${sample}_counts.txt"), emit: stats

    script:

    """
    samtools view -bS ${alignment} > ${sample}.bam
    samtools sort -m ${task.memory.toGiga()}G -@ ${task.cpus} ${sample}.bam -o ${sample}_sorted.bam
    samtools idxstats ${sample}_sorted.bam > ${sample}_counts.txt
    """
}