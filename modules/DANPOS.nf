process danpos_mono{
  container 'uschwartz/danpos'
  memory { params.genomeSize > 200000000 ? params.high_memory : params.low_memory}
  publishDir "${params.outDir}/RUN/01_MONO-NUCS_PROFILE", mode: 'copy', pattern: "*_monoNucs_profile.bw"

  input:
  tuple val(sampleID), file(bam)
  file(chrSizes)

  output:
  file("*_monoNucs_profile.bw")
  tuple val(sampleID), file("result/pooled/*.xls")


  script:
  resolution = ( params.genomeSize > 200000000 ? '10':'1')
  """
  danpos.py dpos $bam -m 1 --extend 70 -c $params.genomeSize \
  -u 0 -z 20 -a $resolution -e 1  > $sampleID"_DANPOS_stats.txt"
  wigToBigWig result/pooled/*.wig -clip $chrSizes $sampleID"_monoNucs_profile.bw"
  """
}

process danpos_sub{
  container 'uschwartz/danpos'
  memory { params.genomeSize > 200000000 ? params.high_memory : params.low_memory}
  publishDir "${params.outDir}/RUN/02_SUB-NUCS_PROFILE", mode: 'copy', pattern: "*_subNucs_profile.bw"

  input:
  tuple val(sampleID), file(bam)
  file(chrSizes)

  output:
  file("*_subNucs_profile.bw")
  tuple val(sampleID), file("result/pooled/*.xls")

  script:
  """
  danpos.py dpos $bam -m 1 --extend 70 -c $params.genomeSize \
  -u 0 -z 70 -a 20 -e 1  > $sampleID"_DANPOS_stats.txt"
  wigToBigWig result/pooled/*.wig -clip $chrSizes $sampleID"_subNucs_profile.bw"
  """
}
