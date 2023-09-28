process statistics_read{
  container 'uschwartz/r_nucmacc:v3.1'
  publishDir "${params.outDir}/QC/06_FragmentStatistics/${sampleID}", mode: 'copy'

  input:
  tuple val(sampleID), file(sieve_mono), file(sieve_sub), file(fastqc), file(alignment), file(qualimap)

  output:
  file("*.txt")
  file("*.pdf")

  script:
  """
  GenerateTxtFragCounts.R
  """
}

process statistics_plot{
  container 'uschwartz/r_nucmacc:v3.1'
  publishDir "${params.outDir}/QC/06_FragmentStatistics", mode: 'copy'

  input:
  file(statistics_read)

  output:
  file("*.txt")
  file("*.pdf")

  script:
  """
  Plot_comparison.R
  """
}
