#!/usr/bin/env nextflow
/*
================================================================================
    annovar-nf
    Scalable, reproducible ANNOVAR annotation of VCF/BCF files (hg19 / hg38)
    https://github.com/assetdaniyarov/annovar-nf
    Author : Asset Daniyarov
    License: MIT

    Modes (--mode):
      annotate  (default) annotate VCF/BCF files
      check_db  read-only audit of the ANNOVAR databases used by the protocols
      update_db install missing/damaged databases and build the newest ClinVar
================================================================================
*/

include { ANNOVAR_NF } from './workflows/annovar_nf'
include { ANNOVAR_DB } from './workflows/annovar_db'
include { helpMessage; pipelineBanner; dbBanner; completionSummary; validateParameters; readInputs; dbBuilds } from './subworkflows/local/utils_annovarnf'

workflow {

    main:
    def ids = []    // input-file ids, for the completion summary
    if (params.help) {
        log.info helpMessage()
    }
    else if (params.mode in ['check_db', 'update_db']) {
        validateParameters()
        log.info dbBanner()
        ANNOVAR_DB(dbBuilds(), params.mode == 'update_db')
    }
    else {
        validateParameters()
        def inputs = readInputs(params.input.toString())
        ids = inputs.collect { e -> e[0].id }
        log.info pipelineBanner(inputs)

        ANNOVAR_NF(channel.fromList(inputs))
    }

    onComplete:
    if (!params.help) {
        log.info completionSummary(ids)
    }
}
