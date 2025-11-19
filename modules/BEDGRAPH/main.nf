process GET_BEGRAPH {
    publishDir "${params.outdir}/Bedgraph/",
        mode: 'copy',
        scratch: true,
        enabled: params.save_tmp
    tag "${bam.baseName}"
    label "lowMemST"
    // conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == "apptainer" ?
        'oras://community.wave.seqera.io/library/deeptools_pygenometracks:225dd30c0028afe4':
        'community.wave.seqera.io/library/deeptools_pygenometracks:c37e66d27ebf153d' }"

    input:
    file(bam)
    file(bai)

    output:
    path("${bam.baseName}.bedgraph"), emit: bedgraph

    script:
    """
    bamCoverage -b ${bam} -o ${bam.baseName}.bedgraph \
        --normalizeUsing CPM --binSize 10 \
        --outFileFormat bedgraph
    wait

    echo "Finished generating bedgraph file for ${bam.baseName}. Output files: ${bam.baseName}.bedgraph"
    """
}




