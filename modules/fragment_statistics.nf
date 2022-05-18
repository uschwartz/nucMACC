process statistics_read{
  container 'uschwartz/r_nucmacc:v3'
  publishDir "${params.outDir}/QC/FragmentStatistics", mode: 'copy'

  input:
  tuple val(sampleID), file(sieve_mono), file(sieve_sub), file(fastqc), file(alignment), file(qualimap)

  output:
  file("*.txt")

  script:
  """
  GenerateTxtFragCounts.R
  """
}

process statistics_plot{
  container 'uschwartz/r_nucmacc:v3'
  publishDir "${params.outDir}/QC/FragmentStatistics", mode: 'copy'

  input:
  

  output:
  file("*.pdf")

  script:
  """
  Plot-comparison.R
  """
}
