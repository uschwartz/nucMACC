process make_TSS_plots_monoNucs{
  container 'uschwartz/r_nucmacc:v3.1'

  publishDir "${params.outDir}/RUN/07_TSS_profile/monoNucs", mode: 'copy'

  input:
  file(input)

  output:
  file("*.pdf")

  script:
  """
  make_TSS_plots_monoNucs.R $input
  """
}

process make_TSS_plots_subNucs{
  container 'uschwartz/r_nucmacc:v3.1'

  publishDir "${params.outDir}/RUN/07_TSS_profile/subNucs", mode: 'copy'

  input:
  file(input)

  output:
  file("*.pdf")

  script:
  """
  make_TSS_plots_subNucs.R $input
  """
}
