process featureCounts_diff_mono{
        label 'big'

        publishDir "${params.outDir}/RUN/02_DIFF_NUCS_READ_COUNTS/mono", mode: 'copy'


        input:
        file(saf)
        file(Nucs)
        val(names)

        output:
        file("diffNucs_readCounts.csv")
        file("diffNucs_readCounts.csv.summary")

        script:
        """
        featureCounts  -F SAF -a $saf \
        -o diffNucs_readCounts.csv \
        --fracOverlap 0.65 \
        -T $task.cpus -p -B -O \
        $Nucs

        names="$names"
        IFS=\$'\t' read -r -a namesArray <<< \$(sed 's/\\[\\|\\]//g; s/,/\t/g' <<< "\$names")
        new_line="Geneid\tChr\tStart\tEnd\tStrand\tLength"
        for name in "\${namesArray[@]}"; do
        new_line+="\t\$name"
        done
        sed -i "2s/.*/\$new_line/" "diffNucs_readCounts.csv"

        """
}

process featureCounts_diff_sub{
        label 'big'

        publishDir "${params.outDir}/RUN/02_DIFF_NUCS_READ_COUNTS/sub", mode: 'copy'


        input:
        file(saf)
        file(Nucs)
        val(names)
        
        output:
        file("diffNucs_readCounts.csv")
        file("diffNucs_readCounts.csv.summary")

        script:
        """
        featureCounts  -F SAF -a $saf \
        -o diffNucs_readCounts.csv \
        --fracOverlap 0.7 \
        -T $task.cpus -p -B -O \
        $Nucs

        names="$names"
        IFS=\$'\t' read -r -a namesArray <<< \$(sed 's/\\[\\|\\]//g; s/,/\t/g' <<< "\$names")
        new_line="Geneid\tChr\tStart\tEnd\tStrand\tLength"
        for name in "\${namesArray[@]}"; do
        new_line+="\t\$name"
        done
        sed -i "2s/.*/\$new_line/" "diffNucs_readCounts.csv"
        """
}
