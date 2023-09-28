process mergeBam_mono{
  label 'mid'
  //publishDir "${params.outDir}/RUN/00_bamEntry/monoNuc", mode: 'copy'

  input:
  tuple val(Sample_Name), path(path_mono)

  output:
  tuple val(Sample_Name), file("*.bam")

  script:
  """
  samtools merge "${Sample_Name}_mono.bam" $path_mono --threads $task.cpus -u
  """
}

process mergeBam_sub{
  label 'mid'
  //publishDir "${params.outDir}/RUN/00_bamEntry/subNuc", mode: 'copy'

  input:
  tuple val(Sample_Name), path(path_sub)

  output:
  tuple val(Sample_Name), file("*.bam")

  script:
  """
  samtools merge "${Sample_Name}_sub.bam" $path_sub --threads $task.cpus -u
  """
}
