include { PYGT_INI as PYGTINI } from '../modules/PYGTINI'
include { TRACKPLOT as TRACKPLOT } from '../modules/pyGenomeTracks'

workflow TRACKPLOT2 {
    take:
    bedgraphs
    peakCallBed
    geneModelGTF
    trackIniNameSuffix
    region_ch
    groups

    main:
    // groups.view()
    // Generate tracks.ini file
    PYGTINI(bedgraphs, peakCallBed, geneModelGTF, trackIniNameSuffix, groups)
    ch_tracks_ini = PYGTINI.out.inifile

    // Generate plots for each tracks.ini file
    def ch_ini_flat = ch_tracks_ini.flatten()
    // ch_ini_flat.view()
    def ch_region_flat = region_ch.flatten()
    // ch_region_flat.view()
    def combined = ch_region_flat.combine(ch_ini_flat)
    // combined.view()
    TRACKPLOT(combined)
    ch_tracks = TRACKPLOT.out.trackplot
    // TRACKPLOT(ch_tracks_ini, region)  //
    // ch_tracks = TRACKPLOT.out.trackplot

    emit:
    tracks_ini_output = ch_tracks_ini
    tracks_output = ch_tracks
}


