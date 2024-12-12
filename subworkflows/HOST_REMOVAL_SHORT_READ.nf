/*
Subworkflow for fastqc from raw data, fastp trimming, then fastqc again
Requries set params:

params.fastp_q  = Q score for trimming

*/


include { BOWTIE2 as BOWTIE2 } from '../modules/BOWTIE2.nf'
include { SAMTOOLS_STATS as SAMTOOLS_STATS } from '../modules/SAMTOOLS_STATS.nf'

workflow HOST_REMOVAL_SHORT_READ {
    take:
        ch_fastqs_trim                // channel: [val(sample), [fastq_1, fastq_2]]
        ch_hostgen                  // channel: [val(sample), fasta]
    main:

        ch_for_bowtie = ch_fastqs_trim.join(ch_hostgen)

        BOWTIE2(ch_for_bowtie)
        SAMTOOLS_STATS(BOWTIE2.out.sam)

}