process alignment{
  label 'big'
  publishDir "${params.outDir}/RUN/02_ALIGNMENT", mode: 'copy', pattern: "*_alignment_stats.txt"


  input:
  tuple val(sampleID), file(read1), file(read2)

  output:
  file "*_alignment_stats.txt"
  file "*_aligned.bam"

  script:
  """
  bowtie2 -t \
  --threads $task.cpus \
  --very-sensitive-local \
  --no-discordant \
  --no-mix \
  --dovetail \
  -x $params.genomeIdx \
  -1 $read1 \
  -2 $read2 \
  2> ${sampleID}_alignment_stats.txt \
  | samtools view -bS -f 2 -@ $task.cpus - > ${sampleID}"_aligned.bam"
  """

}
  // | samtools view -bS -f 2 -@ $task.cpus - > ${sampleID}"_aligned.bam
