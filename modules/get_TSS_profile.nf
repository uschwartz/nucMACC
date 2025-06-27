process TSS_profile_mono{
   if (params.container_engine == 'docker') {
    containerOptions "-v \$(dirname ${params.TSS}):\$(dirname ${params.TSS})"
   }

  label 'big'

  input:
  file(bw)

  output:
  file "computeMatrix2plot_mono.txt.gz"


  script:
  """
  computeMatrix reference-point -S $bw \
   -R $params.TSS \
   --referencePoint TSS \
   -o "computeMatrix2plot_mono.txt.gz" \
   -b 1500 -a 1500 --smartLabels -p $task.cpus


  """
}

process TSS_profile_plot_mono{


  input:
  file(computeMatrix2plot_mono)

  output:
  file "values_Profile_mono.txt"

  script:
  """
  plotProfile -m $computeMatrix2plot_mono \
       -out 'DefaultHeatmap_mono.png' \
       --outFileNameData 'values_Profile_mono.txt'
  """
}

process TSS_profile_sub{
  if (params.container_engine == 'docker') {
    containerOptions "-v \$(dirname ${params.TSS}):\$(dirname ${params.TSS})"
  }
  label 'big'

  input:
  file(bw)

  output:
  file "computeMatrix2plot_sub.txt.gz"


  script:
  """
  computeMatrix reference-point -S $bw \
   -R $params.TSS \
   --referencePoint TSS \
   -o "computeMatrix2plot_sub.txt.gz" \
   -b 1500 -a 1500 --smartLabels -p $task.cpus


  """
}

process TSS_profile_plot_sub{

  input:
  file(computeMatrix2plot_sub)

  output:
  file "values_Profile_sub.txt"

  script:
  """
  plotProfile -m $computeMatrix2plot_sub \
       -out 'DefaultHeatmap_sub.png' \
       --outFileNameData 'values_Profile_sub.txt'

  """
}
