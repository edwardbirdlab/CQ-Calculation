#!/usr/bin/env nextflow

nextflow.enable.dsl=2


if (params.workflow_opt == 'shortread_meta') {

    ch_fastq = Channel.fromPath(params.sample_sheet) \
        | splitCsv(header:true) \
        | map { row-> tuple(row.sample, file(row.r1), file(row.r2)) }

    ch_hostgen = Channel.fromPath(params.sample_sheet) \
        | splitCsv(header:true) \
        | map { row-> tuple(row.sample, file(row.refernce_genome)) }

    }

include { SHORT_READ_METAGENOMIC as SHORT_READ_METAGENOMIC } from './workflows/SHORT_READ_METAGENOMIC.nf'
include { SR_MULTIQC as SR_MULTIQC } from './workflows/SR_MULTIQC.nf'


workflow {


    if (params.workflow_opt == 'shortread_meta') {

        SHORT_READ_METAGENOMIC(ch_fastq, ch_hostgen)

        }

    if (params.workflow_opt == 'sr_multiqc') {

        SR_MULTIQC()

        }

}