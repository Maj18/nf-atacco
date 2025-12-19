process GET_BEGRAPH {
    publishDir "${params.outdir}/topAR/Bedgraph/",
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
    tuple file(bam), file(bai), val(group)

    output:
    tuple path("${bam.baseName}.bedgraph"), val(group), emit: bedgraph

    script:
    """
    bamCoverage -b ${bam} -o ${bam.baseName}.bedgraph \
        --normalizeUsing RPGC --binSize 10 \
        --effectiveGenomeSize 2699716848 \
        --outFileFormat bedgraph \
        --ignoreForNormalization chrX chrY chrY_KI270740v1_random \
        || true

   # Any non-zero exit code becomes 0, so Nextflow won’t fail

    echo "Finished generating bedgraph file for ${bam.baseName}. Output files: ${bam.baseName}.bedgraph"
    """
}




