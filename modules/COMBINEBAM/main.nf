process MERGEBAM {
    publishDir "${params.outdir}/TOBIAS_TFfootprinting/bam_filtered_combined/",
        mode: 'copy'
    tag "${group}"
    label "lowMemMT"
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == "apptainer" ?
        'oras://community.wave.seqera.io/library/bedtools_samtools:384b4ee53b3a9fc6':
        'community.wave.seqera.io/library/bedtools_samtools:73f3ef465b3bd116' }"

    input:
    tuple val(group), file(bams), file(bais)

    output:
    tuple val(group),
          path("${group}_merged.bam"),
          path("${group}_merged.bam.bai"),
          emit: mergedbams_allMAPQ
    tuple val(group),
          path("${group}_merged_filtered.bam"),
          path("${group}_merged_filtered.bam.bai"),
          emit: mergedbams

    script:
    """
    echo "Merging bams files for ${group}..."
    samtools merge -@ "${task.cpus}" "${group}_merged.bam" ${bams.join(' ')}
    sync
    samtools index -@ "${task.cpus}" "${group}_merged.bam"

    # Only keep reads with MAPQ>30
    samtools view -q 30 -b -f 2 "${group}_merged.bam" > "${group}_merged_filtered.bam"
    sync

    echo "Creating index for ${group}_merged.bam..."
    samtools index -@ "${task.cpus}" "${group}_merged_filtered.bam"

    echo "Processing of ${group} completed. Merged bam and index files generated."
    """
}
