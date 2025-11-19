process PYGT_INI {
    publishDir "${params.outdir}/TrackVisualization/",
        mode: 'copy',
        scratch: true,
        enabled: params.save_tmp
    tag "track .ini"
    label "lowMemST"

    input:
    val(bedgraphs)
    val(peakCallBed)
    val(geneModelGTF)
    val(trackIniNameSuffix)

    output:
    path("tracks${trackIniNameSuffix}.ini"), emit: inifile

    script:
    """
    bash ${moduleDir}/templates/writeini.sh ${bedgraphs} ${peakCallBed} ${geneModelGTF} ${trackIniNameSuffix}
    """
}


