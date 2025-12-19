process ATACORRECT {
    publishDir "${params.outdir}/TOBIAS_TFfootprinting/ATACorrect/${group}/",
        mode: 'copy'
    tag "${group}"
    label "highMemMT"
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == "apptainer" ?
        'oras://community.wave.seqera.io/library/tobias:0.17.2--f96a578ea350926b':
        'community.wave.seqera.io/library/tobias:0.17.2--3a49d50188327fd2' }"

    input:
    tuple val(group), val(group_peak), file(merged_bam), file(merged_bai)

    output:
    tuple val(group), path("corrected"), emit: corrected_dir

    script:
    """
    echo "Before running ATACorrect, please make sure your bam files have been filtered (e.g. MAPQ>30)," 
    echo "please also make sure your peak files have been filtered as well to only keep peaks with coverage > 100 reads & qvalue<0.01!"
    echo "Run ATACorrect for ${group}..."
    echo "${group_peak}"
    head "${group_peak}"   
    TOBIAS ATACorrect \
        --bam "${merged_bam}" \
        --genome "${params.refgenome}" \
        --peaks "${group_peak}" \
        --outdir "corrected" \
        --cores "${task.cpus}"

    echo "Correcting ${group} completed."
    """
}

