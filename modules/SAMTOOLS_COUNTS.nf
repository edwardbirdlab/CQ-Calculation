process SAMTOOLS_COUNTS {
    label 'samtoolssort'
    container 'ebird013/samtools:1.17'

    input:
        tuple val(sample), file(alignment)
    output:
        tuple val(sample), path("${sample}_counts.txt"), emit: stats

    script:

    """
    samtools sort -n -m ${task.memory.toGiga()}G -@ ${task.cpus} ${alignment} -o ${sample}_sorted.sam
    samtools idxstats ${sample}_sorted.sam > ${sample}_counts.txt
    """
}