process sieve_mono{
  label 'big'
  publishDir "${params.outDir}/RUN/04_ALIGNMENT_FILTERING/monoNuc", mode: 'copy', pattern: "*_mono_FiltLog.txt"

  input:
  tuple val(sampleID), val(bam)

  output:
  tuple val(sampleID), file("*_mono_FiltLog.txt")
  tuple val(sampleID), file("*_mono.bam")

  script:
  blacklistOpt = ( params.blacklist ? "--blackListFileName $params.blacklist":'')
  """
  alignmentSieve -b $bam \
  -o ${sampleID}"_mono.bam" \
  -p $task.cpus \
  --filterMetrics  ${sampleID}"_mono_FiltLog.txt" \
  --minFragmentLength 140 \
  --maxFragmentLength 200 \
  $blacklistOpt
  """
}

process sieve_sub{
  label 'big'
  publishDir "${params.outDir}/RUN/04_ALIGNMENT_FILTERING/subNuc", mode: 'copy', pattern: "*_sub_FiltLog.txt"

  input:
  tuple val(sampleID), val(bam)

  output:
  tuple val(sampleID), file("*_sub_FiltLog.txt")
  tuple val(sampleID), file("*_sub.bam")

  script:
  blacklistOpt = ( params.blacklist ? "--blackListFileName $params.blacklist":'')
  """
  alignmentSieve -b $bam \
  -o ${sampleID}"_sub.bam" \
  -p $task.cpus \
  --filterMetrics  ${sampleID}"_sub_FiltLog.txt" \
  --maxFragmentLength 139 \
  $blacklistOpt
  """
}
