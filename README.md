# nucMACC

**Original paper** - https://doi.org/10.1101/2022.12.29.521985

## Introduction
nucMACC is an automated analysis pipeline for the analysis of nucleosome positions, accessibility and stability. The pipeline contains two main workflows:

1. `MNaseQC` for QC and  exploratory analysis
2. `nucMACC` for analysis of nucleosome positions, accessibility and stability

Given trimmed paired-end sequencing reads in fastq format, this pipeline will run:

* `MNaseQC` and `nucMACC`
    1. QC using `FastQC` on fastq files
    2. Alignment using `Bowtie2` on fastq files
    3. QC using `Qualimap` on aligned fragments
    4. Fragment size distribution plot
    5. Group the fragments by size using 'deepTools alignmentSieve' and optionally filter blacklisted regions
        1. Mono-nucleosome (140 - 200 bp)
        2. Sub-nucleosome (< 140 bp)
    6. Report fragment statistics of each processing step
    7. Create nucleosome maps of Mono- and Sub-nucleosomes using `DANPOS`
    8. Optionally create TSS profiles using `deepTools`
    9. Summary reports using `MultiQC`


* `MNaseQC` specific
    1. PCA of nucleosome maps using `deepTools`
    2. Correlation analysis using `deepTools`

* `nucMACC` specific
    1. Pool all mono-nucleosome samples
    2. Obtain mono-nucleosome positions from pooled samples and sub-nucleosome positions from lowest MNase concentration using `DANPOS`
    3. Get GC content of nucleosome positions using `bedtools genomecov`
    4. Remove nucleosome positions with low fragment count (mono-nucleosmes < 30 and sub-nucleosomes < 5)
    5. Calculate (sub-)nucMACC score using linear regression Analysis
    6. Correct for MNase GC-bias using LOWESS
    7. Identify hyper-/hypo-accessible nucleosomes or unstable and non-canoncical nucleosomes.

`nucMACC` is meant to run on pooled replicates in fastq format, whereas `MNaseQC` uses single replicates. As the `MNaseQC` and the `nucMACC` workflow have several steps in common, it is recommended to run first `MNaseQC` and report the fragment size selected bam files using `--publishBamFlt`. Then setting `--bamEntry` option, a shorter version of the `nucMACC` workflow can be run using the generated bam files as input. Here in an additional step at the beginning replicates are pooled.

## Get started

### Requirements

`Docker` and `nextflow` are required to run the nucMACC pipeline. Additional software used in the pipeline is packaged in Docker container and will be automatically downloaded during the first execution of the pipeline.
The pipeline is compatible with all computational infrastructures. Executing the pipeline on cloud or HPC systems may require to adapt the [`nextflow.config`](https://www.nextflow.io/docs/latest/basic.html).
For large reference genomes the pipeline requires at least 32 GB memory and default settings allocate 45 GB memory to demanding processes. Memory usage can be adjusted using the option `--high_memory` or directly in the [`nextflow.config`](https://www.nextflow.io/docs/latest/basic.html).

### Installation
You can obtain the pipeline directly from GitHub:
```bash
git clone https://github.com/uschwartz/nucMACC.git
```


## Contact

Please log all issues/suggestions on the nucMACC GitHub page: https://github.com/uschwartz/nucMACC/issues

Uwe Schwartz: uwe.schwartz@ur.de

## Cite

Wernig-Zorc et al., 2023, nucMACC: An optimized MNase-seq pipeline measures genome-wide nucleosome accessibility and stability. bioRxiv (https://doi.org/10.1101/2022.12.29.521985)
