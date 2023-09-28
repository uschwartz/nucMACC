process qualimap {

  container 'uschwartz/qualimap'
  label 'mid'
  publishDir "${params.outDir}/QC/03_QUALIMAP", mode: 'copy'

  input:
  tuple val(sampleID), file(bam)

  output:
  file "${sampleID}"
  tuple val(sampleID), file("${sampleID}/*.txt")


  script:
  """
  qualimap bamqc --java-mem-size=14G -bam $bam -c -outdir ${sampleID}
  """
}
