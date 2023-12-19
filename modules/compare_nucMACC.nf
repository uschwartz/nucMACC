process compare_nucMACC_mono{
container 'leoschmutterer/diff_nucmacc:v1'
publishDir "${params.outDir}/RUN/04_DIFF_NUCS_COMP/", mode: 'copy', pattern: "*.tsv"
publishDir "${params.outDir}/RUN/04_DIFF_NUCS_COMP/REGIONS_MONO", mode: 'copy', pattern: "*.bed"

input:
file(input)
file(diff_bed)
file(diff_stats)

output:
file("*.tsv")
file("*.bed")
script:
"""
compare_nucMACC_mono.R $input $diff_bed $diff_stats
"""
}

process compare_nucMACC_sub{
container 'leoschmutterer/diff_nucmacc:v1'
publishDir "${params.outDir}/RUN/04_DIFF_NUCS_COMP/", mode: 'copy', pattern: "*.tsv"
publishDir "${params.outDir}/RUN/04_DIFF_NUCS_COMP/REGIONS_SUB", mode: 'copy', pattern: "*.bed"
input:
file(input)
file(diff_bed)
file(diff_stats)

output:
file("*.tsv")
file("*.bed")
script:
"""
compare_nucMACC_sub.R $input $diff_bed $diff_stats
"""
}