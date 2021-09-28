process pool{
  label 'big'

  input:
  file(monos)
  output:
  tuple val("pooled"), file("pooled_monoNucs.bam")
  file("chrom_Sizes.txt")

  script:
  """
  samtools merge -@ 6 pooled_monoNucs.bam $monos
  samtools view -H pooled_monoNucs.bam \
  | awk -v OFS='\t' '/^@SQ/ {print \$2,\$3}' \
  | sed  's/.N://g' >chrom_Sizes.txt
  """
}
