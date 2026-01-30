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

// include { FILTERBAM } from '../modules/FILTERBAM'
include { GET_BEGRAPH2  } from './subworkflows/BEDGRAPH.nf'
include { TRACKPLOT2 } from './subworkflows/TRACK_VISUALIZATION.nf'
include { TOBIAS } from './subworkflows/FOOTPRINTING_TOBIAS.nf'
include { TFmotifEnrichment } from './subworkflows/MONALISA_TFMOTIFENRICHMENT.nf'
include { DECOUPLER_DIFFTFACTIVITY } from './subworkflows/DECOUPLER_DIFFTFACTIVITY.nf'
include { TF_INTEGRATION } from './subworkflows/TF_INTEGRATION.nf'

if (params.sampleSheet) {
    ch_input = Channel
        .fromPath(params.sampleSheet)
        .splitCsv(header: true, sep: ',')
        .map { row -> tuple(file(row.file),file(row.bai),row.group) }
}

workflow ENTRY_TRACKPLOT {
    // Plot genome browser track plots for provide peak regions:
    if ( !params.peakCallBedfile || !params.geneModelGTFfile || !params.regions ) {
        error "Please provide --peakCallBedfile --geneModelGTFfile --regions for TRACKPLOT analysis"
    }
    // Prepare the bedgraph files
    if( params.bedgraphFiles ) {
       // User provides bedgraph files
       // Sort the bedgraph files based on group
       ch_sorted = Channel.fromPath(params.bedgraphFiles)
            .splitCsv(header: true, sep: ',')
            .map { row -> tuple(file(row.file),row.group) }
	    .groupTuple(by: 1)
	    .transpose()
       ch_sorted.view()
       bedgraphs_input = ch_sorted
            .map { it[0].toString() }
            .collect()
            .map { it.join(',')}
       group_input = ch_sorted
            .map { it[1].toString() }
            .collect()
            .map { it.join(',') }
       // bedgraphs_input.view()
    } else {
        // Generate bedgraphs from BAMs via GET_BEGRAPH2
        GET_BEGRAPH2(ch_input)
        // Sort the bedgraph files based on groups
        ch_sorted = GET_BEGRAPH2.out.bedgraph_output
	        .groupTuple(by: 1)
            .transpose()
        // ch_sorted.view()
        // Collect into one list channel
        bedgraphs_input = ch_sorted
            .map { it[0].toString() }      // convert each file to path string
            .collect()                  // gather all emitted paths into a List
            .map { it.join(',') }       // join the list with commas
        bedgraphs_input.view()
        group_input = ch_sorted
            .map { it[1].toString() }      // convert each file to path string
            .collect()                  // gather all emitted paths into a List
            .map { it.join(',') }       // join the list with commas
	group_input.view()
    }

    // Make genome browser track plots for each provide regions
    // sanity check / required error
    if( !params.bedgraphFiles && params.nobg ) {
        error "You must provide --bedgraphFiles when --nobg is set"
    }
    peakCallBed_input = Channel.fromPath(params.peakCallBedfile)
    geneModelGTF_input = Channel.fromPath(params.geneModelGTFfile)
    trackIniNameSuffix_input = Channel.value(params.trackIniNameSuffix)
    region_input = Channel.from(params.regions.split(','))
    // region_input.view()
    TRACKPLOT2(bedgraphs_input, 
        peakCallBed_input, 
        geneModelGTF_input, 
        trackIniNameSuffix_input,
        region_input,
        group_input)
}


workflow ENTRY_TOBIAS {
    // Run TOBIAS TF footprinting analysis
    if ( !params.sampleSheet || !params.group_peaks || !params.peakAnnotation || !params.refgenome ) {
        error "Please provide --sampleSheet --group_peaks --peakAnnotation --refgenome for TOBIAS footprinting analysis"
    }
    // Filtering the bam files to only keep reads with MAPQ>30
    // FILTERBAM(file_input)
    // ch_filteredbam = FILTERBAM.out.filteredbam
    // ch_groupBams = ch_filteredbam.map { bam, bai, group -> tuple(group, bam, bai) }
    // The bam fiels will be filtered after merging replicates in TOBIAS module
    ch_groupBams = ch_input.map { bam, bai, group -> tuple(group, bam, bai) }
        .groupTuple(by: 0)                 // group by group name
    ch_groupPeaks = Channel.from(params.group_peaks.split(','))
        .map { it.split(':') }
        .map { parts -> tuple(parts[0], parts[1]) }
    ch_peakAnnotation = Channel.fromPath(params.peakAnnotation)
    ch_groupBamsPeaks = ch_groupBams
        .join(ch_groupPeaks)
        .combine(ch_peakAnnotation)
        // .map { joined, peakAnnotation -> tuple(joined[0], joined[1], joined[2], joined[3], peakAnnotation) }.view() //
    TOBIAS(ch_groupBamsPeaks)
}


workflow ENTRY_MONALISA {
    // Run MONALISA TF motif enrichment analysis
    if ( !params.difftables || !params.peakAnnotation ) {
        error "Please provide both --difftables and --peakAnnotation for MONALISA_TFmotifEnrichment"
    }
    ch_difftable = Channel.value(params.difftables)
    ch_peakAnnotation = Channel.fromPath(params.peakAnnotation)
    ch_peakAnnotation.view()
    if ( !params.pfm_file ) {
        ch_pfm_file = Channel.fromPath("${projectDir}/data/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt")
    } else {
        ch_pfm_file = Channel.fromPath(params.pfm_file)
    }
    ch_pfm_file.view()
    TFmotifEnrichment(ch_difftable, ch_peakAnnotation, ch_pfm_file)
}

workflow ENTRY_DIFFTFACTIVITY {
    def has_dds = params.dds
    def has_log = params.logNormCount && params.design
    if ( !has_dds && !has_log ) {
        error "Please provide either --dds OR both --logNormCount and --design, for ENTRY_DIFFTFACTIVITY"
    }
    if ( has_dds && has_log ) {
        error "Please provide EITHER --dds OR --logNormCount + --design, not both, for ENTRY_DIFFTFACTIVITY"
    }
    if ( !params.group_order ) {
        error "Please also provide --group_order for ENTRY_DIFFTFACTIVITY"
    }

    // Run decoupleR differential TF activity analaysis
    DECOUPLER_DIFFTFACTIVITY()
}

workflow ENTRY_TFINTEGRATION {
    // Run TF integration
    ch_diffTFexpr_files = Channel.value(params.diffTFexpr_files)
    ch_diffTFbinding_file = Channel.value(params.diffTFbinding_file)
    ch_diffTFactivity_files = Channel.value(params.diffTFactivity_files)
    ch_monalisa_files = Channel.value(params.monalisa_files)
    TF_INTEGRATION(ch_diffTFexpr_files,
        ch_diffTFbinding_file, ch_diffTFactivity_files, ch_monalisa_files)
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

