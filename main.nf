#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/atacco
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : 
    Contact: Yuan.li@nbis.se
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

println """\
         ATACCO   P I P E L I N E
         ===================================
         GitHub: 
         ___________________________________
         SAMPLE SHEET   : ${params.sampleSheet}
         OUTPUT DIR     : ${params.outdir}
        ___________________________________
         """
         .stripIndent()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GET_BEGRAPH2  } from './subworkflows/BEDGRAPH.nf'
include { TRACKPLOT2 } from './subworkflows/TRACK_VISUALIZATION.nf'

if (params.sampleSheet) {
    ch_input = Channel
        .fromPath(params.sampleSheet)
        .splitCsv(header: true, sep: ',')
        .map { row -> tuple(file(row.file),file(row.bai),row.group) }
} else if (params.file) {
    ch_input = Channel
        .fromPath(params.file)
        .map { f -> tuple(file(f), 'default_group') }
} else {
    error "Please provide either --sampleSheet or --file"
}

workflow {
    if( params.bedgraphFiles ) {
        // Case 1: user provides bedgraph files
       println "bedgraphFiles param: ${params.bedgraphFiles}"
       // bedgraphs_input = Channel.from(params.bedgraphFiles.split(',')).map { it.toString() }.collect().map { it.join(' ') }
       bedgraphs_input = Channel.from(params.bedgraphFiles)
       bedgraphs_input.view()
    } else {
        // Case 2: generate bedgraphs from BAMs via GET_BEGRAPH2
        file_input = ch_input.map { it[0] }
        bai_input  = ch_input.map { it[1] }

        def bg = GET_BEGRAPH2(file_input, bai_input)
        // Collect into one list channel
        // To do, this may need to be tested
        bedgraphs_input = bg.bedgraph_output
                            .map { it.toString() }      // convert each file to path string
                            .collect()                  // gather all emitted paths into a List
                            .map { it.join(',') }       // join the list with commas
        bedgraphs_input.view()
    }

    // sanity check / required error
    if( !params.bedgraphFiles && params.nobg ) {
        error "You must provide --bedgraphFiles when --nobg is set"
    }
    peakCallBed_input = Channel.fromPath(params.peakCallBedfile)
    geneModelGTF_input = Channel.fromPath(params.geneModelGTFfile)
    trackIniNameSuffix_input = Channel.value(params.trackIniNameSuffix)
    region_input = Channel.from(params.regions.split(','))
    // region_input.view()
    TRACKPLOT2(bedgraphs_input, peakCallBed_input, geneModelGTF_input, trackIniNameSuffix_input, region_input)
}

workflow.onComplete {
    println( workflow.success ? """
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """ : """
        Failed: ${workflow.errorReport}
        exit status : ${workflow.exitStatus}
        """
    )
}

