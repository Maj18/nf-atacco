process DIFFTFACTIVITY {
    publishDir "${params.outdir}/",
        mode: 'copy'
    tag "Dorothea_VIPER"
    label "lowMemST"
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == "apptainer" ?
        'oras://community.wave.seqera.io/library/bioconductor-biocparallel_bioconductor-bsgenome.hsapiens.ucsc.hg38_bioconductor-complexheatmap_bioconductor-deseq2_pruned:db3861e19b7a781c':
        'community.wave.seqera.io/library/bioconductor-biocparallel_bioconductor-bsgenome.hsapiens.ucsc.hg38_bioconductor-complexheatmap_bioconductor-deseq2_pruned:75a508058f405a9a' }"

    output:
    path("diffTFactivity_RNA"), emit: diffTFactivity_dir

    script:
    """
    Rscript ${moduleDir}/templates/decoupleR.R --dds "${params.dds}" --outdir "diffTFactivity_RNA" --group_order "${params.group_order}" --currentCovariate "${params.currentCovariate}"
    """
}

