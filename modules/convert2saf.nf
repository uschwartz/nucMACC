process convert2saf_mono{
  publishDir "${params.outDir}/RUN/07_MONO-NUCS_POSITIONS", mode: 'copy', pattern: "*_nucPositions.bed"

  input:
  tuple val(sampleID), file(xls), file(bam)

  output:
  file("*_nucPositions.bed")

  script:
  """
  echo -e "NucID\tChr\tStart\tEnd\tStrand\tPeak_score" > $sampleID"_nucPositions.saf"
  awk -v OFS='\t' 'NR > 1  {print "nuc"NR,\$1,\$2,\$3,".",\$5}' \
   $xls >> $sampleID"_nucPositions.saf"
  awk -v OFS='\t' 'NR > 1 {print \$1,\$2,\$3,"nuc"NR,\$5,"."}' $xls \
  >> $sampleID"_nucPositions.bed"
  """
}
