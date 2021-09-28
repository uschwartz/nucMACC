process danpos_mono{
  container 'uschwartz/danpos'
  input:
  tuple val(sampleID), val(bam)

  //output:

  script:
  """
  danpos.py dpos $bam -m 1 --extend 70 -c 1000000 \
   -u 0 -z 1 -a 1 -e 1  > $sampleID"_DANPOS_stats.txt"
  """
}
