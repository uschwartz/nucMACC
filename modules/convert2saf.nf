process convert2saf_mono{
  publishDir "${params.outDir}/RUN/07_NUCS_POSITIONS", mode: 'copy', pattern: "*_nucPositions.bed"

  input:
  tuple val(sampleID), file(xls), file(bam)

  output:
  file("*_nucPositions.bed")
  file("*_nucPositions.saf")

  script:
  """
  echo -e "GeneID\tChr\tStart\tEnd\tStrand\tPeak_score" > $sampleID"_nucPositions.saf"
  awk -v OFS='\t' 'NR > 1  {print "nuc"NR,\$1,\$2,\$3,".",\$5}' \
   $xls >> $sampleID"_nucPositions.saf"
  awk -v OFS='\t' 'NR > 1 {print \$1,\$2,\$3,"nuc"NR,\$5,"."}' $xls \
  >> $sampleID"_nucPositions.bed"
  """
}

process convert2saf_sub{
  publishDir "${params.outDir}/RUN/07_NUCS_POSITIONS", mode: 'copy', pattern: "*_sub-nucPositions.bed"

  input:
  tuple val(sampleID), file(xls), file(bam)

  output:
  file("*_sub-nucPositions.bed")
  file("*_sub-nucPositions.saf")

  script:
  """
  echo -e "GeneID\tChr\tStart\tEnd\tStrand\tPeak_score" > $sampleID"_sub-nucPositions.saf"
  awk -v OFS='\t' 'NR > 1  {print "nuc"NR,\$1,\$2,\$3,".",\$5}' \
   $xls >> $sampleID"_sub-nucPositions.saf"
  awk -v OFS='\t' 'NR > 1 {print \$1,\$2,\$3,"nuc"NR,\$5,"."}' $xls \
  >> $sampleID"_sub-nucPositions.bed"
  """
}
