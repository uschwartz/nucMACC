process plot_correlation_Heatmap{
  publishDir "${params.outDir}/QC/07_correlationHeatmap", mode: 'copy'

  input:
  file(bigwigSummary)
  val(name)

  output:
  file "*.pdf"

  script:
  """
  plotCorrelation -in $bigwigSummary \
    --corMethod $params.correlationMethod --skipZeros \
    --plotTitle "${params.correlationMethod} correlation ${name}" \
    --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
    --plotFile "${name}_correlationHeatmap_${params.correlationMethod}.pdf"   \
  """
}
