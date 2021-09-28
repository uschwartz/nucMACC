process pool{
  label 'big'

  input:
  file(monos)
  output:
  tuple val("pooled"), file("pooled_monoNucs.bam")

  script:
  """
  samtools merge -@ 6 pooled_monoNucs.bam $monos
  """
}
