process MONALISA {
    publishDir "${params.outdir}/",
        mode: 'copy'
    tag "Monalisa"
    label "lowMemST"
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == "apptainer" ?
        'oras://community.wave.seqera.io/library/bioconductor-biocparallel_bioconductor-bsgenome.hsapiens.ucsc.hg38_bioconductor-complexheatmap_bioconductor-genomicranges_pruned:9e8dadbd18083cb4':
        'community.wave.seqera.io/library/bioconductor-biocparallel_bioconductor-bsgenome.hsapiens.ucsc.hg38_bioconductor-complexheatmap_bioconductor-genomicranges_pruned:d96030d6d371c514' }"

    input:
    val(difftables)
    file(peakAnnotation)
    file(pfm_file)

    output:
    path("MonaLisa"), emit: MonaLisaBinnedMotifEnrichment

    script:
    """
    Rscript ${moduleDir}/templates/monalisaBinEnr.R "${difftables}" ${peakAnnotation} ${pfm_file} "MonaLisa/" "${task.cpus}"
    """
}

