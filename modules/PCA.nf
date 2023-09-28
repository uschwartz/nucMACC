process plot_PCA{
  publishDir "${params.outDir}/QC/08_PCA", mode: 'copy'

  input:
  file(bigwigSummary)
  val(name)

  output:
  file "${name}_PCA.pdf"

  script:
  """
  plotPCA -in $bigwigSummary \
  --plotFile ${name}_PCA.pdf \
  --plotTitle "PCA ${name}" \
  --transpose \
  """
}
