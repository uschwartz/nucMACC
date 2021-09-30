process featureCounts_mono{
        label 'big'

        publishDir "${params.outDir}/RUN/09_Nuc_readCounts", mode: 'copy'


        input:
        file(saf)
        file(monoNucs)

        output:
        file("monoNucs_readCounts.csv")

        script:
        """
        featureCounts  -F SAF -a $saf \
        -o monoNucs_readCounts.csv \
        -T $task.cpus -p -B --largestOverlap \
        $monoNucs
        """
}
