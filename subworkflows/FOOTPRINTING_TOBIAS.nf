include { MERGEBAM } from '../modules/COMBINEBAM'
include { ATACORRECT } from '../modules/TOBIAS_ATACORRECT'
include { FOOTPRINTSCORES } from '../modules/TOBIAS_FOOTPRINTSCORES'
include { BINDETECT } from '../modules/TOBIAS_BINDETECT'
include { FOOTPRINTPLOT } from '../modules/TOBIAS_FOOTPRINTPLOT'
include { BINDINGHEATMAP } from '../modules/TOBIAS_TFBINDINGHEATMAP'

workflow TOBIAS {
    take:
    bamspeaks_grouped

    main:
    // Merge bams by groups
    bams_grouped = 
        bamspeaks_grouped.map { group, bams, bais, peaks, peakAnnotation -> tuple(group, bams, bais) }.view() //
    MERGEBAM(bams_grouped)
    ch_mergedbam = MERGEBAM.out.mergedbams

    // TOBIAS ATACorrect
    // This will do Tn5 shifting by default
    // It will also correct Tn5 binding sequence preferences
    peaks_grouped = bamspeaks_grouped.map { group, bams, bais, peaks, peakAnnotation -> tuple(group, peaks) } //
    megedbam_peaks = peaks_grouped.join(ch_mergedbam).view()
    ATACORRECT(megedbam_peaks)
    ch_corrected_dir = ATACORRECT.out.corrected_dir

    // TOBIAS FootprintFootprintScores
    correctedSignal_peaks = peaks_grouped.join(ch_corrected_dir).view()
    FOOTPRINTSCORES(correctedSignal_peaks)
    ch_ftscores = FOOTPRINTSCORES.out.ftscores
    // ch_ftscores.view()

    // TOBIAS BINDetect
    groups = ch_ftscores.map { it[0].toString() }      // convert each file to path string
                            .collect()                  // gather all emitted paths into a List
                            .map { it.join(',') }       // join the list with commas
    groups.view()
    ftscores = ch_ftscores.map { it[1].toString() }      // convert each file to path string
                            .collect()                  // gather all emitted paths into a List
                            .map { it.join(',') }       // join the list with commas
    peakAnn = bamspeaks_grouped
        .map { group, bams, bais, peaks, peakAnnotation -> peakAnnotation }
        .unique()
    // ftscores.view()
    // groups_ftscores = groups.combine(ftscores).combine(peakAnn)
    groups_ftscores = groups
        .combine(ftscores)
        .combine(peakAnn)
        // .map { grpFtscore, peakAnnotation -> tuple(grpFtscore[0], grpFtscore[1], peakAnnotation) }
    groups_ftscores.view()
    BINDETECT(groups_ftscores)
    ch_diffTFbinding = BINDETECT.out.DiffTFBinding_dir

    // TOBIAS footprint plot
    TFs = ch_diffTFbinding
        .map { folder ->
            folder.listFiles()
                .findAll { it.isDirectory() }
                .collect { it.name }
        }
        .flatten()
    // TFs = Channel.of(params.TFs)
    //     .map { it.split(",") }
    //     .flatten()
    TFs.view()
    FTplotIn_ch = TFs.combine(ch_diffTFbinding).combine(groups).combine(
        ch_corrected_dir.map { it[1].toString() } 
            .collect()              // gather all emitted paths into a List
            .map { it.join(',') })
    FTplotIn_ch.view()
    FOOTPRINTPLOT(FTplotIn_ch)
    ch_footprint_plot = FOOTPRINTPLOT.out.footprint_plot

    // TOBIAS TF binding heatmap
    BINDINGHEATMAP(FTplotIn_ch)
    ch_binding_heatmap = BINDINGHEATMAP.out.binding_heatmap

    emit:
    ch_corrected_dir = ch_corrected_dir
    ch_ftscores = ch_ftscores
    ch_diffTFbinding = ch_diffTFbinding
    ch_footprint_plot = ch_footprint_plot
    ch_binding_heatmap = ch_binding_heatmap
}

