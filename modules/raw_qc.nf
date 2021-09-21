process fastqc{

  publishDir "${params.outDir}/RUN/01_FASTQC/${sampleID}", mode: 'copy', pattern: "*.html"

  input:
  tuple val(sampleID), file(read)

  output:
  file "*_fastqc.zip"
  file "*_fastqc.html"

  script:
  """
  fastqc $read
  """
}

process multiqc{

  publishDir "${params.outDir}/QC/multiqc/", mode: 'copy'

  input:
  file('*')
  output:
  file "*multiqc_report.html"
  file "*_data"

  script:
  """
  multiqc .
  """
}
