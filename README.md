# nucMACC

**Original paper** - https://doi.org/10.1101/2022.12.29.521985

## Introduction
nucMACC is an automated analysis pipeline for the analysis of nucleosome positions, accessibility and stability. The pipeline contains two main workflows:

1. `MNaseQC` for QC and  exploratory analysis
2. `nucMACC` for analysis of nucleosome positions, accessibility and stability

Given trimmed paired-end sequencing reads in fastq format, this pipeline will run:

* `MNaseQC` and `nucMACC`
        * QC using `FastQC` on fastq files
        * Alignment using `Bowtie2` on fastq files
        * QC using `Qualimap` on aligned fragments
        * Fragment size distribution plot
        * Group the fragments by size using 'deepTools alignmentSieve' and optionally filter blacklisted regions
                * Mono-nucleosome (140 - 200 bp)
                * Sub-nucleosome (< 140 bp)
        * Report fragment statistics of each processing step          
        * Create nucleosome maps of Mono- and Sub-nucleosomes using `DANPOS`
        * Optionally create TSS profiles using `deepTools`
        * Summary reports using `MultiQC`

* `MNaseQC` specific
        * PCA of nucleosome maps using `deepTools`
        * Correlation analysis using `deepTools`

* `nucMACC` specific
        * Pool all mono-nucleosome samples
        * Obtain mono-nucleosome positions from pooled samples and sub-nucleosome positions from lowest MNase concentration using `DANPOS`
        * Get GC content of nucleosome positions using `bedtools genomecov`
        * Remove nucleosome positions with low fragment count (mono-nucleosmes < 30 and sub-nucleosomes < 5)
        * Calculate (sub-)nucMACC score using linear regression Analysis
        * Correct for MNase GC-bias using LOWESS
        * Identify hyper-/hypo-accessible nucleosomes or unstable and non-canoncical nucleosomes.

`nucMACC` is meant to run on pooled replicates in fastq format, whereas `MNaseQC` uses single replicates. As the `MNaseQC` and the `nucMACC` workflow have several steps in common it is recommended to run first `MNaseQC` and report the grouped bam files using `--publishBamFlt`. Then setting `--bamEntry` option, a shorter version of the `nucMACC` workflow can be run using the generated bam files as input. Here in an additional step at the beginning replicates are pooled.

## Contact

Please log all issues/suggestions on the nucMACC GitHub page: https://github.com/uschwartz/nucMACC/issues

Uwe Schwartz: uwe.schwartz@ur.de

## Cite

Wernig-Zorc et al., 2023, nucMACC: An optimized MNase-seq pipeline measures genome-wide nucleosome accessibility and stability. bioRxiv (https://doi.org/10.1101/2022.12.29.521985)
