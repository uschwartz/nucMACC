process diff_nucMACC_mono{
    container 'leoschmutterer/diff_nucmacc:v1.1'
    publishDir "${params.outDir}/RUN/03_DIFF_NUCS_POS/mono", mode: 'copy'

    input:
    file(csv)
    file(counts)
    file(summary)

    output:
    file("Diff_nucMACC_result_table.tsv")
    file("Diff_nucMACC_fdrFilt_result_table.tsv")
    file("Positions_Diff_nucMACC_fdrFilt.bed")
    file("*_scored_positions.bed")

    script:
    """
    diff_nucMACC_EdgeR.R $csv $counts $summary
    """

}


process diff_nucMACC_sub{
    container 'leoschmutterer/diff_nucmacc:v1.1'
    publishDir "${params.outDir}/RUN/03_DIFF_NUCS_POS/sub", mode: 'copy'
    
    input:
    file(csv)
    file(counts)
    file(summary)
    
    output:
    file("Diff_nucMACC_result_table.tsv")
    file("Diff_nucMACC_fdrFilt_result_table.tsv")
    file("Positions_Diff_nucMACC_fdrFilt.bed")
    file("*_scored_positions.bed")
    script:
    """
    diff_nucMACC_EdgeR.R $csv $counts $summary
    """
}
