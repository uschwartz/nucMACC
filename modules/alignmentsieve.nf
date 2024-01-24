process sieve_mono{
  if(params.blacklist){
  containerOptions "-v \$(dirname ${params.blacklist}):\$(dirname ${params.blacklist})"
  }
  label 'big'
  publishDir "${params.outDir}/QC/05_ALIGNMENT_FILTERING/monoNuc", mode: 'copy', pattern: "*_mono_FiltLog.txt"
  publishDir "${params.outDir}/RUN/00_ALIGNMENT/monoNuc", mode: 'copy', pattern: "*_mono.bam", enabled:params.publishBamFlt

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
  if(params.blacklist){
  containerOptions "-v \$(dirname ${params.blacklist}):\$(dirname ${params.blacklist})"
  }
  label 'big'
  publishDir "${params.outDir}/QC/05_ALIGNMENT_FILTERING/subNuc", mode: 'copy', pattern: "*_sub_FiltLog.txt"
  publishDir "${params.outDir}/RUN/00_ALIGNMENT/subNuc", mode: 'copy', pattern: "*_sub.bam", enabled:params.publishBamFlt

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
