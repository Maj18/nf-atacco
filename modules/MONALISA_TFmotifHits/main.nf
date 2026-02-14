process TFMOTIFHITS {
    publishDir "${params.outdir}/",
        mode: 'copy'
    tag "MotifHits"
    label "lowMemST"
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == "apptainer" ?
        'oras://community.wave.seqera.io/library/bioconductor-biocparallel_bioconductor-bsgenome.hsapiens.ucsc.hg38_bioconductor-complexheatmap_bioconductor-genomicranges_pruned:9e8dadbd18083cb4':
        'community.wave.seqera.io/library/bioconductor-biocparallel_bioconductor-bsgenome.hsapiens.ucsc.hg38_bioconductor-complexheatmap_bioconductor-genomicranges_pruned:d96030d6d371c514' }"

    input:
        val(difftables)
        file(peakAnnotation)
        val(selectedPeaks)
        file(pfm_file)

    output:
        path("PeakTFmotifHits/TFmotif_hits_list_allPeaks.RDS"), emit: hits_list_rds
        path("PeakTFmotifHits/TFmotif_hitsMatrix_allPeaks.RDS"), emit: hits_matrix_rds
        path("PeakTFmotifHits/TFmotif_hitsMatrix_allPeaks.tsv"), emit: hits_matrix_tsv
        path("PeakTFmotifHits/TFmotif_hitsMatrix_diff/"), emit: hits_matrix_diff_tsv_dir
        path("PeakTFmotifHits/Monalisa_sessionInfo.txt"), emit: sessionInfo

    script:
    """
    Rscript ${moduleDir}/templates/monalisaTFmotifHits.R "${difftables}" ${peakAnnotation} "${selectedPeaks}" ${pfm_file} "${task.cpus}"
    """
}

