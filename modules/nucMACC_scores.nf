process nucMACC_scores{
  container 'uschwartz/r_nucmacc:v3.1'

  publishDir "${params.outDir}/RUN/05_nucMACC", mode: 'copy'

  input:
  file(readCounts)
  file(input)
  file(counts_sum_mono)

  output:
  file("Figures/*.png")
  tuple file("hyperAcc_monoNucs.bed"), file("hypoAcc_monoNucs.bed")
  file("nucMACC_result_table.tsv")
  file("nucMACC_scores.bedgraph")
  file("positions_monoNucs.bed")


  script:
  """
  get_nucMACC_scores.R $input $readCounts $counts_sum_mono
  """
}

process sub_nucMACC_scores{
  container 'uschwartz/r_nucmacc:v3.1'

  publishDir "${params.outDir}/RUN/06_sub-nucMACC", mode: 'copy'

  input:
  file(readCounts)
  file(input)
  val(minConc)
  file(nucMACC_table)
  file(readCounts_nucs)
  file(counts_sum_mono)
  file(counts_sum_sub)

  output:
  file("Figures/*.png")
  tuple file("nonCanonical_subNucs.bed"), file("unStable_subNucs.bed")
  file("sub-nucMACC_result_table.tsv")
  file("sub-nucMACC_scores.bedgraph")
  file("positions_subNucs.bed")


  script:
  """
  get_sub-nucMACC_scores.R $input $readCounts $minConc $nucMACC_table $readCounts_nucs $counts_sum_mono $counts_sum_sub
  """
}
