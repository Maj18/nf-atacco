process TRACKPLOT {
    publishDir "${params.outdir}/TrackVisualization/",
        mode: 'copy',
        scratch: true,
        enabled: params.save_tmp
    tag "pyGenomeeracks"
    label "lowMemST"
    // conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == "apptainer" ?
        'oras://community.wave.seqera.io/library/deeptools_pygenometracks:225dd30c0028afe4':
        'community.wave.seqera.io/library/deeptools_pygenometracks:c37e66d27ebf153d' }"

    input:
    tuple val(region), file(ini)
    // file(ini)
    // val(region) // e.g. chr14:91295621-95725599
    // def (region2, name) = $region.split(';')

    output:
    path("*.png"), emit: trackplot

    script:
    // strip suffix and leading "tracks"
    def (region2, name) = region.split(';')
    def name2 = name.replace(':', '_').replace('-', '_').replace('(', '_').replace(')', '')
    println "region = $region2"
    println "name = $name2"
    """
    pyGenomeTracks --tracks ${ini} --region ${region2} --trackLabelFraction 0.2 --dpi 130 -o ${name2}.png

    echo "Finished generating track plot for ${name2}. Output files: ${name2}.png"
    """
}

