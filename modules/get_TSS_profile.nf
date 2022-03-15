process TSS_profile_mono{

  publishDir "${params.outDir}/RUN/11_TSS_profile/monoNucs", mode: 'copy'

  input:
  file(bw)

  output:
  file "computeMatrix2plot_mono.txt.gz"
  file "computeMatrix2txt_mono.txt.gz"
  file "computeMatrix_geneList_mono.bed"

  script:
  """
  computeMatrix reference-point -S $bw \
   -R $params.TSS \
   --referencePoint TSS \
   -o "computeMatrix2plot_mono.txt.gz" \
   --outFileNameMatrix "computeMatrix2txt_mono.txt.gz" \
   --outFileSortedRegions "computeMatrix_geneList_mono.bed" \
   -b 1500 -a 1500 --smartLabels -p 10


  """
}

process TSS_profile_plot_mono{

  publishDir "${params.outDir}/RUN/11_TSS_profile/monoNucs", mode: 'copy'

  input:
  file(computeMatrix2plot_mono)

  output:
  file "*.png"
  file "sortedRegions_Profile_mono.txt"
  file "values_Profile_mono.txt"

  script:
  """
  plotProfile -m $computeMatrix2plot_mono \
       -out 'DefaultHeatmap_mono.png' \
       --outFileSortedRegions 'sortedRegions_Profile_mono.txt'\
       --outFileNameData 'values_Profile_mono.txt'

  plotProfile -m $computeMatrix2plot_mono \
        -out 'DefaultHeatmap_mono_grouped.png' \
        --perGroup

  """
}

process TSS_profile_sub{

  publishDir "${params.outDir}/RUN/11_TSS_profile/subNucs", mode: 'copy'

  input:
  file(bw)

  output:
  file "computeMatrix2plot_sub.txt.gz"
  file "computeMatrix2txt_sub.txt.gz"
  file "computeMatrix_geneList_sub.bed"

  script:
  """
  computeMatrix reference-point -S $bw \
   -R $params.TSS \
   --referencePoint TSS \
   -o "computeMatrix2plot_sub.txt.gz" \
   --outFileNameMatrix "computeMatrix2txt_sub.txt.gz" \
   --outFileSortedRegions "computeMatrix_geneList_sub.bed" \
   -b 1500 -a 1500 --smartLabels -p 10


  """
}

process TSS_profile_plot_sub{

  publishDir "${params.outDir}/RUN/11_TSS_profile/subNucs", mode: 'copy'

  input:
  file(computeMatrix2plot_sub)

  output:
  file "*.png"
  file "sortedRegions_Profile_sub.txt"
  file "values_Profile_sub.txt"

  script:
  """
  plotProfile -m $computeMatrix2plot_sub \
       -out 'DefaultHeatmap_sub.png' \
       --outFileSortedRegions 'sortedRegions_Profile_sub.txt'\
       --outFileNameData 'values_Profile_sub.txt'

  plotProfile -m $computeMatrix2plot_sub \
        -out 'DefaultHeatmap_sub_grouped.png' \
        --perGroup

  """
}
