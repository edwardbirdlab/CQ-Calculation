/*
Subworkflow for fastqc from raw data, fastp trimming, then fastqc again
Requries set params:

params.fastp_q  = Q score for trimming

*/


include { BOWTIE2 as BOWTIE2 } from '../modules/BOWTIE2.nf'
include { SAMTOOLS_STATS as SAMTOOLS_STATS } from '../modules/SAMTOOLS_STATS.nf'
include { SAMTOOLS_COUNTS as SAMTOOLS_COUNTS } from '../modules/SAMTOOLS_COUNTS.nf'
include { SAMTOOLS_SAM2BAM as SAMTOOLS_SAM2BAM } from '../modules/SAMTOOLS_SAM2BAM.nf'
include { SAMTOOLS_BAM_SORT_POS as SAMTOOLS_BAM_SORT_POS } from '../modules/SAMTOOLS_BAM_SORT_POS.nf'


workflow HOST_REMOVAL_SHORT_READ {
    take:
        ch_fastqs_trim                // channel: [val(sample), [fastq_1, fastq_2]]
        ch_hostgen                  // channel: [val(sample), fasta]
    main:

        ch_for_bowtie = ch_fastqs_trim.join(ch_hostgen)

        BOWTIE2(ch_for_bowtie)
        SAMTOOLS_STATS(BOWTIE2.out.sam)
        SAMTOOLS_SAM2BAM(BOWTIE2.out.sam)
        SAMTOOLS_BAM_SORT_POS(SAMTOOLS_SAM2BAM.out.bam)
        SAMTOOLS_COUNTS(SAMTOOLS_BAM_SORT_POS.out.sort)

}