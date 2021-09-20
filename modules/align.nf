process alignment{
  label 'big'
  echo true

  input:
  tuple val(sampleID), file(read1), file(read2)

  //output:

  script:
  """
  bowtie2 --version
  """

}
