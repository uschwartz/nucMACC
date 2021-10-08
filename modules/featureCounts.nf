process featureCounts_mono{
        label 'big'

        publishDir "${params.outDir}/RUN/08_NUCS_READ_COUNTS", mode: 'copy'


        input:
        file(saf)
        file(monoNucs)

        output:
        file("monoNucs_readCounts_wGC.csv")

        script:
        """
        featureCounts  -F SAF -a $saf \
        -o monoNucs_readCounts.csv \
        -T $task.cpus -p -B --largestOverlap \
        $monoNucs

        #get number of columns
        numb=\$(cat monoNucs_readCounts.csv | awk 'NR > 2 {print NF; exit}')
        last_c=\$((\$numb-1))
        gc_c=\$((\$numb+1))


        awk -v OFS='\t' 'NR > 2 {for(i=2;i<=NF;++i)printf \$i""FS; print}' monoNucs_readCounts.csv | cut -f 2- \
        | bedtools nuc -fi $params.genome -bed - | cut -f 1-\$last_c,\$gc_c \
        > monoNucs_readCounts_wGC.csv


        """
}


process featureCounts_sub{
        label 'big'

        publishDir "${params.outDir}/RUN/8_NUCS_READ_COUNTS", mode: 'copy'


        input:
        file(saf)
        file(subNucs)

        output:
        file("subNucs_readCounts_wGC.csv")

        script:
        """
        featureCounts  -F SAF -a $saf \
        -o subNucs_readCounts.csv \
        -T $task.cpus -p -B --largestOverlap \
        $subNucs

        #get number of columns
        numb=\$(cat subNucs_readCounts.csv | awk 'NR > 2 {print NF; exit}')
        last_c=\$((\$numb-1))
        gc_c=\$((\$numb+1))


        awk -v OFS='\t' 'NR > 2 {for(i=2;i<=NF;++i)printf \$i""FS; print}' subNucs_readCounts.csv | cut -f 2- \
        | bedtools nuc -fi $params.genome -bed - | cut -f 1-\$last_c,\$gc_c \
        > subNucs_readCounts_wGC.csv


        """
}
