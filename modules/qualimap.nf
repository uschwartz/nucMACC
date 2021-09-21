process qualimap {
  echo true
  //label 'big'
  //publishDir "${params.outDir}/RUN/02_ALIGNMENT", mode: 'copy', pattern: "*_alignment_stats.txt"

  input:
  tuple val(sampleID), file(bam)

  //output:

  script:
  """
  qualimap --help
  """
}

//qualimap bamqc --java-mem-size=16G -bam $bam -c -outdir ${sampleID} 
