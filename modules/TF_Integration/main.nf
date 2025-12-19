process INTEGRATION {
    publishDir "${params.outdir}/",
        mode: 'copy'
    tag "integration"
    label "lowMemST"
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == "apptainer" ?
        'oras://community.wave.seqera.io/library/bioconductor-biocparallel_bioconductor-bsgenome.hsapiens.ucsc.hg38_bioconductor-complexheatmap_bioconductor-deseq2_pruned:a1ceaff5bb99ce9d':
        'community.wave.seqera.io/library/bioconductor-biocparallel_bioconductor-bsgenome.hsapiens.ucsc.hg38_bioconductor-complexheatmap_bioconductor-deseq2_pruned:7a508c273d8b955a' }"

    input:
    val(diffTFexpr_files)
    val(diffTFbinding_file)
    val(diffTFactivity_files)
    val(monalisa_files)

    output:
    path("Integration"), emit: Integration_dir

    script:
    """
    Rscript ${moduleDir}/templates/Integration.R --outdir "Integration" --diffTFexpr_files "$diffTFexpr_files" --diffTFbinding_file "$diffTFbinding_file" --diffTFactivity_files "$diffTFactivity_files" --monalisa_files "$monalisa_files" --group_order "${params.group_order}" --human_mouse_symbol_ortholog "${projectDir}/data/human_mouse_symbol_ortholog.csv" --ggvenn_script "${moduleDir}/templates/ggvenn.R"
    """
}


