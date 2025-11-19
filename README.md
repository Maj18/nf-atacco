# nf-atacco

| ![Workflow](https://img.shields.io/badge/Workflow-Nextflow-blue) | ![Version](https://img.shields.io/badge/Version-1.0.0-green) |
|------------------------------------------------------------------|-------------------------------------------------------------|

### Pipeline Stats
- **Modules**: 3 (BEDGRAPH, PYGTINI, pyGenomeTracks)
- **Languages**: Nextflow DSL2, Groovy, R
- **Supported Platforms**: Local, HPC, Cloud
- **Execution Profiles**: Local, Docker, Apptainer, UPPMAX

### Quick Links
- [Documentation](https://github.com/Maj18/nf-atacco#readme)
- [Issues](https://github.com/Maj18/nf-atacco/issues)
- [Nextflow](https://www.nextflow.io/)
- [NBIS Epigenomics Workshop](https://nbis-workshop-epigenomics.readthedocs.io/)
- [Seqera](https://seqera.io/)

`nf-atacco` is a bioinformatics pipeline built using [Nextflow](https://www.nextflow.io/) for ATAC-seq (Assay for Transposase-Accessible Chromatin using sequencing) downstream analysis. It currently covers bam to bedgraph conversion, and genome track visualizaiton.

*`nf-atacco` is currently in development and more modules will be added and optimizations performed in time.*


## Features

- **Reproducibility**: Ensures consistent results across runs using Nextflow.
- **Scalability**: Supports execution on local machines, HPC clusters, and cloud platforms.
- **Customizable**: Easily configurable to suit specific experimental needs.

## Requirements

- [Nextflow](https://www.nextflow.io/)
- [Docker](https://www.docker.com/) or [Apptainer](https://sylabs.io/) (optional for containerized execution)

## Installation

1. Clone the repository:
    ```bash
    git clone https://github.com/Maj18/nf-atacco.git
    cd nf-atacco
    ```

2. Install nf-core
    ```bash
    https://nf-co.re/docs/nf-core-tools/installation
    ```

3. Ensure dependencies are installed (e.g., Docker/Singularity).

### Offline setup

If you need to run the pipeline in an offline environment, you can download all the container images to the `conts` directory. You can find the container information in the module/*/main.nf file, use singularity pull ... or docker pull ... to download the containers.

After this, run the pipeline with the `offline` profile in addition to your preferred profile(s).


### HPC profiles (SUPR NAISS Sweden only)

`pdc_kth` profile for PDC Dardel cluster
`uppmax` profile for UPPMAX cluster including auto-detect for `miarka`, `pelle`, `snowy` and `rackham` queues.

If you are running the pipeline on these HPC clusters, you can use the specific profile. These profiles are configured to use Singularity containers and has specific settings for the HPC environments.
Profiles originally written by `Pontus Freyhult (@pontus)` and adapted for the pipeline by `Aditya Singh (@addityea)`.

## Usage

Run the pipeline with the following command:

```bash
nextflow run main.nf --project <project id> \
    --sampleSheet <sampleSheet> --outdir <outDir> -profile uppmax,offline \
    --peakCallBedfile <peakCallBedfile_merged> \
    --geneModelGTFfile <geneModelGTFfile> \
    --regions <regions> \
    --save_tmp true \
    -with-report report.html -with-trace trace.txt \
    -with-timeline timeline.html -resume -bg  

```

Example regions 'chr1:156591000-156592814;chr1:genrich_9654(NAXE),chr3:51705518-53705890;chr3:genrich_99114(SPCS1)'

You may customise the pipeline by modifying the `nextflow.config` file or by passing additional parameters in the command line.
There are some profiles available for different environments:
- 'local' for running on a local machine
- 'docker' for running with Docker
- 'arm' for specific flags for ARM architecture
- 'apptainer' for running with Apptainer
- 'UPPMAX' for running on UPPMAX cluster

## Custom Entry points

The pipeline provides several custom entry points to execute specific parts of the workflow without running the entire pipeline. To be added.

To use a specific entry point, run the pipeline with the `-entry` flag followed by the desired entry point name. For example:


## Parameters

The pipeline accepts the following parameters, which can be customized to suit your analysis needs:

- **`--sampleSheet`**: Path to the sample sheet file. *(Default: `null`)*
- **`--outdir`**: Directory where output files will be saved. *(Default: `results`)*
- **`--cpus`**: Number of threads to use for computation. *(Default: `14`)*
- **`--project`**: UPPMAX account to use for job submission. *(Default: `null`)*
- **`--memory`**: Amount of memory allocated per process. *(Default: `20 GB`)*
- **`--time`**: Maximum runtime for processes. *(Default: `10-00:00:00`)*
- **`--retry`**: Number of retries for failed processes. *(Default: `3`)*
- **`--highMemForks`**: Maximum number of high-memory forks. *(Default: `2`)*
- **`--fatMemForks`**: Maximum number of fat-memory forks. *(Default: `1`)*

Example of a sample sheet:

```csv
file,group
/path/to/BAM/sample_A.bam,Control
/path/to/BAM/sample_B.bam,Control
/path/to/BAM/sample_C.bam,Treatment
/path/to/BAM/sample_D.bam,Treatment
```

#### Optional Flags
- **`--nobg`**: skip convert bam to bedgraph. *(Default: `false`)*

#### Convert filtered bam (should have corresponding .bai file available took) to bedgraph
- **`--bedgraphFiles`**: .bedgraph files, e.g. (A1.bedgraph A2.bedgraph) *(Default: `null`)*
- **`--peakCallBedfile`**: consensus peak .bed file *(Default: `null`)*
- **`--geneModelGTFfile`**: reference gtf file *(Default: `null`)*

#### Peak Calling
- **`--trackIniNameSuffix`**: name suffix for the tracks.ini file needed by pyGenomeTracks *(Default: `null`)*
- **`--regions`**: genome regions to be plotted s*(Default: `chr1:1000000-2000000,chr2:1000000-2000000`)*


## Output

The pipeline generates the following outputs:
- bedgraph files
- genome browser track plots


## Programs used

| Program | Purpose | Link |
| ------- | ------- | ---- |
| deepTools | processing, normalizing, comparing, and visualizing NGS data | [deepTools](https://deeptools.readthedocs.io/) |
| pyGenomeTracks | produce high-quality genome browser tracks | [pyGenomeTracks](https://pygenometracks.readthedocs.io/en/latest/) |


## Contributing

Contributions are welcome! Please fork the repository and submit a pull request.

## Acknowledgments

- Built with [Nextflow](https://www.nextflow.io/)
- Inspired by the [NBIS Epigenomics Workshop repo](https://nbis-workshop-epigenomics.readthedocs.io/en/latest/content/tutorials/atacseq_tutorials.html)
- Used [addityea/nf-atac](https://github.com/addityea/nf-atac.git) as a template



