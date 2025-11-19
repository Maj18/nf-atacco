include { GET_BEGRAPH as BEGRAPH } from '../modules/BEDGRAPH'

workflow GET_BEGRAPH2 {
    take:
    file_input
    bai_input

    main:
    BEGRAPH(file_input, bai_input)
    ch_bedgraph = BEGRAPH.out.bedgraph 

    emit:
    begraph_output = ch_bedgraph
}
