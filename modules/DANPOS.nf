process danpos_mono{
  container 'uschwartz/danpos'

  publishDir "${params.outDir}/RUN/05_MONO-NUCS_PROFILE", mode: 'copy', pattern: "*_monoNucs_profile.bw"

  input:
  tuple val(sampleID), val(bam)
  file(chrSizes)

  output:
  file("*_monoNucs_profile.bw")

  script:
  """
  danpos.py dpos $bam -m 1 --extend 70 -c $params.genomeSize \
  -u 0 -z 1 -a 1 -e 1  > $sampleID"_DANPOS_stats.txt"
  wigToBigWig result/pooled/*.wig -clip $chrSizes $sampleID"_monoNucs_profile.bw"
  """
}

process danpos_sub{
  container 'uschwartz/danpos'

  publishDir "${params.outDir}/RUN/06_SUB-NUCS_PROFILE", mode: 'copy', pattern: "*_subNucs_profile.bw"

  input:
  tuple val(sampleID), val(bam)
  file(chrSizes)

  output:
  file("*_subNucs_profile.bw")

  script:
  """
  danpos.py dpos $bam -m 1 --extend 70 -c $params.genomeSize \
  -u 0 -z 70 -a 20 -e 1  > $sampleID"_DANPOS_stats.txt"
  wigToBigWig result/pooled/*.wig -clip $chrSizes $sampleID"_subNucs_profile.bw"
  """
}
