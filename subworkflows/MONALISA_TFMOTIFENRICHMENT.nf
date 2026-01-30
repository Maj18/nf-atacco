include { MONALISA } from '../modules/MONALISA_TFmotifEnrichment'

workflow TFmotifEnrichment {
    take:
    difftable
    peakAnnotation
    pfm_file

    main:
    // Binned motif enrichment
    MONALISA(difftable, peakAnnotation, pfm_file)
    ch_hits_list_rds = MONALISA.out.hits_list_rds
    ch_hits_matrix_rds = MONALISA.out.hits_matrix_rds
    ch_hits_matrix_tsv = MONALISA.out.hits_matrix_tsv
    ch_hits_matrix_diff_tsv = MONALISA.out.hits_matrix_diff_tsv
    ch_binnedmotifenr = MONALISA.out.MonaLisaBinnedMotifEnrichment

    emit:
    ch_hits_list_rds
    ch_hits_matrix_rds
    ch_hits_matrix_tsv
    ch_hits_matrix_diff_tsv
    ch_binnedmotifenr
}

