process InsertSize_Histogram{
  container 'uschwartz/r_nucmacc:v3.1'

  publishDir "${params.outDir}/QC/04_FragmentSize_profile", mode: 'copy'

  input:
  file('*')

  output:
  file("*.pdf")

  script:
  """
  InsertSize_Histogram.R
  """
}
