process FOOTPRINTSCORES {
    publishDir "${params.outdir}/TOBIAS_TFfootprinting/FootprintScores/${group}/",
        mode: 'copy'
    tag "${group}"
    label "highMemMT"
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == "apptainer" ?
        'oras://community.wave.seqera.io/library/tobias:0.17.2--f96a578ea350926b':
        'community.wave.seqera.io/library/tobias:0.17.2--3a49d50188327fd2' }"

    input:
    tuple val(group), val(group_peak), path(dir)

    output:
    tuple val(group), file("${group}_footprint.bw"), emit: ftscores

    script:
    """
    ls -la ${dir}
    echo "Calculating footprint scores for ${group}..."
    TOBIAS FootprintScores \
        --signal "${dir}/${group}_merged_filtered_corrected.bw" \
        --regions "${group_peak}" \
        --output "${group}_footprint.bw" \
        --cores "${task.cpus}"

    echo "Footprint score calculation for ${group} completed."
    """
}

