include { INTEGRATION } from '../modules/TF_Integration'

workflow TF_INTEGRATION {
    take:
    diffTFexpr_files
    diffTFbinding_file
    diffTFactivity_files
    monalisa_files

    main:
    INTEGRATION(diffTFexpr_files, diffTFbinding_file, diffTFactivity_files, monalisa_files)
    ch_integration = INTEGRATION.out.Integration_dir

    emit:
    ch_integration = ch_integration
}
