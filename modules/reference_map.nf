process reference_map_mono{
  container 'uschwartz/r_core:v4.2'
  publishDir "${params.outDir}/RUN/01_REFERENCE_MAP/mono", mode: 'copy'

  input:
  file(input)

  output:
  file ("NucPosRef_allNucs.bed")
  file ("*.saf")

  script:
  """
  NucReferenceMaps_mono.R $input

  awk 'OFS="\t" {print \$1"."\$2"."\$3, \$1, \$2, \$3, "."}' NucPosRef_allNucs.bed > NucPosRef_allNucs.saf
  """
}

process reference_map_sub{
  container 'uschwartz/r_nucmacc:v3.1'
  publishDir "${params.outDir}/RUN/01_REFERENCE_MAP/sub", mode: 'copy'

  input:
  file(input)

  output:
  file ("NucPosRef_allNucs.bed")
  file ("*.saf")

  script:
  """
  NucReferenceMaps_sub.R $input

  awk 'OFS="\t" {print \$1"."\$2"."\$3, \$1, \$2, \$3, "."}' NucPosRef_allNucs.bed > NucPosRef_allNucs.saf
  """
}
