include { MONALISA } from '../modules/MONALISA_TFmotifEnrichment'

workflow TFmotifEnrichment {
    take:
    difftable
    peakAnnotation
    pfm_file

    main:
    // Binned motif enrichment
    MONALISA(difftable, peakAnnotation, pfm_file)
    ch_binnedmotifenr = MONALISA.out.MonaLisaBinnedMotifEnrichment

    emit:
    ch_binnedmotifenr = ch_binnedmotifenr
}

