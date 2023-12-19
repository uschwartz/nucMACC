process bigwigSummary{
  label 'mid'

  input:
  file(bw)
  val(name)

  output:
  //data for correlationHeatmap and PCA
  file "${name}_bw_summary.npz"
  //raw data
  file "${name}_rawCounts.tab"

  script:
  regionOpt = (params.test ? "--region chr3R:18750000:23750000":'')
  """
  multiBigwigSummary bins -b $bw \
  $regionOpt \
  --smartLabels \
  --outFileName ${name}"_bw_summary.npz" \
  --outRawCounts ${name}"_rawCounts.tab" \
  -p $task.cpus
  """
}
