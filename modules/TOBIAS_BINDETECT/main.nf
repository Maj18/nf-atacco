process BINDETECT {
    publishDir "${params.outdir}/TOBIAS_TFfootprinting/BINDetect/",
        mode: 'copy'
    tag "bindetect"
    label "highMemMT"
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == "apptainer" ?
        'oras://community.wave.seqera.io/library/tobias:0.17.2--f96a578ea350926b':
        'community.wave.seqera.io/library/tobias:0.17.2--3a49d50188327fd2' }"

    input:
    tuple val(groups), val(ftscores)
    
    output:
    path("DiffTFBinding"), emit: DiffTFBinding_dir

    script:
    """
    # Extract header for the annotated peaks
    cut -f 1-16,18 "${params.group_peak_annotated}" | head -n 1 > "AnnotatedPeaks_all_header.txt"

    # Match the header and the annotated peaks
    cut -f 1-16,18 "${params.group_peak_annotated}" > "AnnotatedPeaks_all2.bed"

    echo "Filter for HUMAN (!!!!) standard chromosomes only (removes contigs and alternative assemblies)"
    grep -E '^chr[0-22XY]+[[:space:]]' "AnnotatedPeaks_all2.bed" > "AnnotatedPeaks_all3.bed"
    echo ${groups.split(',').join(" ")}
    echo "Identify differential transcription factor (TF) binding between conditions..."
    TOBIAS BINDetect \
        --motifs "${projectDir}/data/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt" \
        --signals ${ftscores.split(',').join(" ")} \
        --genome "${params.refgenome}" \
        --peaks "AnnotatedPeaks_all3.bed" \
        --peak_header "AnnotatedPeaks_all_header.txt" \
        --outdir "DiffTFBinding" \
        --cond_names ${groups.split(',').join(" ")} \
        --cores "${task.cpus}"

    echo "Differential TF binding analysis completed."
    """
}
