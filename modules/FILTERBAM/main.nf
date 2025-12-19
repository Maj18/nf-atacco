process FILTERBAM {
    publishDir "${params.outdir}/topAR/bam_filtered_MAPQ30/",
        mode: 'copy',
        scratch: true,
        enabled: params.save_tmp
    tag "${bam.baseName}"
    label "lowMemST"
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == "apptainer" ?
        'oras://community.wave.seqera.io/library/bedtools_samtools:384b4ee53b3a9fc6':
        'community.wave.seqera.io/library/bedtools_samtools:73f3ef465b3bd116' }"

    input:
    tuple file(bam), file(bai), val(group)

    output:
    tuple file("${bam.baseName}_filtered_MAPQ30.bam"), 
        file("${bam.baseName}_filtered_MAPQ30.bam.bai"), 
        val(group), 
        emit: filteredbam

    script:
    """
    echo "Only keep reads with MAPQ > 30"
    samtools view -q 30 -b -f 2 "${bam}" > "${bam.baseName}_filtered_MAPQ30.bam"
    sync

    echo "Index the filtered bam files"
    samtools index -@ "${task.cpus}" "${bam.baseName}_filtered_MAPQ30.bam"

    echo "Finished filtering and indexing for ${bam.baseName}."
    """
}

