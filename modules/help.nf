/*
* Help message
*/

def helpMessage() {
    println ''
    log.info """
    nucMACC   P I P E L I N E
    =============================
    Usage:

    nextflow run uschwartz/nucMACC --csvInput 'path2csvFile' --outDir 'path2outDir' --genomeIdx 'path2bowtie2_idx' --genomeSize 'eff. genome size' --genome 'path2ref_genome_fasta'

    Mandatory arguments:
      --csvInput        [string] Path to comma-separated file containing information about the samples in the experiment (see ./toyData/input.csv as example) template provided in ./input_template.csv
      --genomeIdx       [string] Path and prefix of bowtie2 index (minus .X.bt2)
      --genomeSize      [integer] Effective genome size, defined as the length of the mappable genome. Used for normalisation (default: 162367812 (dm3))
      --genome          [string] Path to reference genome in fasta Format

    optional arguments:
      --outDir          [string] Name of output directory, which will be created (default: ~/nucMACC_test/)
      --blacklist       [string] A BED file containing regions that should be excluded from all nucleosome analysis (default: false)
      --TSS             [string] A Transcript annotation file in GTF format to obtain TSS profiles (default: false)
      --publishBam      [boolean]if set, aligned bam files will be written to outDir (default: false)
      --publishBamFlt   [boolean]if set, size selected bam files will be written to outDir (default: false)

     """.stripIndent()
     println ''
}
