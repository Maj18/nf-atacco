process BINDINGHEATMAP {
    publishDir "${params.outdir}/TOBIAS_TFfootprinting/BINDetect/DiffTFBinding/${TF}/plots/",
        mode: 'copy'
    tag "heatmap"
    label "lowMemST"
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == "apptainer" ?
        'oras://community.wave.seqera.io/library/tobias:0.17.2--f96a578ea350926b':
        'community.wave.seqera.io/library/tobias:0.17.2--3a49d50188327fd2' }"

    input:
    tuple val(TF), path("DiffTFBinding"), val(groups), val(corrected_dirs)


    output:
    path("${TF}_bindingHeatmap_*vs*.png"), emit: binding_heatmap

    script:
    """
    bash ${moduleDir}/templates/TFbindingheatmap.sh ${TF} ${groups} ${corrected_dirs}
    """
}
