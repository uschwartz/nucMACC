process alignment{
  label 'bigCPU'
  memory { params.genomeSize > 200000000 ? params.high_memory : params.low_memory}
  publishDir "${params.outDir}/QC/02_ALIGNMENT", mode: 'copy', pattern: "*_alignment_stats.txt"
  publishDir "${params.outDir}/RUN/00_ALIGNMENT", mode: 'copy', pattern: "*_aligned.bam", enabled:params.publishBam


  input:
  tuple val(sampleID), file(read1), file(read2)

  output:
  file "*_alignment_stats.txt"
  tuple val(sampleID), file("*_aligned.bam")
  tuple val(sampleID), file("*_alignment_stats.txt")

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
  | samtools view -bS -q 30 -f 2 -@ $task.cpus - | samtools sort -@ $task.cpus - > ${sampleID}"_aligned.bam"

  samtools index -b ${sampleID}"_aligned.bam"
  """

}
