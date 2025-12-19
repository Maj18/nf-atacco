include { GET_BEGRAPH as BEGRAPH } from '../modules/BEDGRAPH'
include { FILTERBAM } from '../modules/FILTERBAM'

workflow GET_BEGRAPH2 {
    take:
    file_input

    main:
    // Filtering the bam files to only keep reads with MAPQ>30
    FILTERBAM(file_input)
    ch_filteredbam = FILTERBAM.out.filteredbam
    // Get Bedgraph files
    BEGRAPH(ch_filteredbam)
    ch_bedgraph = BEGRAPH.out.bedgraph 

    emit:
    bedgraph_output = ch_bedgraph
}
