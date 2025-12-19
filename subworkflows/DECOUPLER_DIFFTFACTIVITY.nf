include { DIFFTFACTIVITY } from '../modules/DECOUPLER_diffTFactivity_RNA'

workflow DECOUPLER_DIFFTFACTIVITY {
    main:
    DIFFTFACTIVITY()
    ch_difftfactivity = DIFFTFACTIVITY.out.diffTFactivity_dir

    emit:
    ch_difftfactivity = ch_difftfactivity
}
