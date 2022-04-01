process InsertSize_Histogram{
  container 'uschwartz/r_nucmacc:v2'

  publishDir "${params.outDir}/QC/FragmentSize_profile", mode: 'copy'

  input:
  file('*')

  output:
  file("*.pdf")

  script:
  """
  InsertSize_Histogram.R
  """
}
