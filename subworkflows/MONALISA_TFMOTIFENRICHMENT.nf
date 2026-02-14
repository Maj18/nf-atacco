include { MONALISA } from '../modules/MONALISA_TFmotifEnrichment'
include { TFMOTIFHITS } from '../modules/MONALISA_TFmotifHits'

workflow TFmotifEnrichment {
    take:
    difftable
    peakAnnotation
    selectedPeaks
    pfm_file

    main:
    // Binned motif enrichment
    MONALISA(difftable, peakAnnotation, pfm_file)
    ch_binnedmotifenr = MONALISA.out.MonaLisaBinnedMotifEnrichment

    TFMOTIFHITS(difftable, peakAnnotation, selectedPeaks, pfm_file)
    ch_hits_list_rds = TFMOTIFHITS.out.hits_list_rds
    ch_hits_matrix_rds = TFMOTIFHITS.out.hits_matrix_rds
    ch_hits_matrix_tsv = TFMOTIFHITS.out.hits_matrix_tsv
    ch_hits_matrix_diff_tsv_dir = TFMOTIFHITS.out.hits_matrix_diff_tsv_dir
    sessionInfo = TFMOTIFHITS.out.sessionInfo

    emit:
    ch_binnedmotifenr
    ch_hits_list_rds
    ch_hits_matrix_rds
    ch_hits_matrix_tsv
    ch_hits_matrix_diff_tsv_dir
    sessionInfo
}

