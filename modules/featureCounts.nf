process featureCounts_mono{
        label 'big'

        publishDir "${params.outDir}/RUN/08_NUCS_READ_COUNTS", mode: 'copy'


        input:
        file(saf)
        file(monoNucs)

        output:
        file("monoNucs_readCounts_wGC.csv")
        file("monoNucs_readCounts.csv.summary")

        script:
        """
        featureCounts  -F SAF -a $saf \
        -o monoNucs_readCounts.csv \
        --fracOverlap 0.7 \
        -T $task.cpus -p -B --largestOverlap \
        $monoNucs

        #get number of columns
        numb=\$(cat monoNucs_readCounts.csv | awk 'NR > 2 {print NF; exit}')
        last_c=\$((\$numb-1))
        gc_c=\$((\$numb+1))


        #get GC content of called nucleosomes
        awk -v OFS='\t' 'NR > 2 {for(i=2;i<=NF;++i)printf \$i""FS; print}' monoNucs_readCounts.csv | cut -f 2- \
        | bedtools nuc -fi $params.genome -bed - | cut -f 1-\$last_c,\$gc_c \
        | awk -v OFS='\t' 'NR > 1 {print}'  > pre_monoNucs_readCounts_wGC.csv

        #add nucID
        awk -v OFS='\t' 'NR > 2 {print \$1}' monoNucs_readCounts.csv \
        | paste pre_monoNucs_readCounts_wGC.csv - >pre_monoNucs_readCounts_wGC_ext.csv

        #prepare header
        awk -v OFS='\t' 'FNR == 2 {print}' monoNucs_readCounts.csv | cut -f 2-\$numb > header_featureCounts.csv
        echo "GC_cont\tnucID" | paste header_featureCounts.csv - > header_all.csv

        #add header
        cat header_all.csv  pre_monoNucs_readCounts_wGC_ext.csv > monoNucs_readCounts_wGC.csv

        """
}


process featureCounts_sub{
        label 'big'

        publishDir "${params.outDir}/RUN/08_NUCS_READ_COUNTS", mode: 'copy'


        input:
        file(saf)
        file(subNucs)

        output:
        file("subNucs_readCounts_wGC.csv")
        file("subNucs_readCounts.csv.summary")

        script:
        """
        featureCounts  -F SAF -a $saf \
        -o subNucs_readCounts.csv \
        --fracOverlap 0.7 \
        -T $task.cpus -p -B --largestOverlap \
        $subNucs

        #get number of columns
        numb=\$(cat subNucs_readCounts.csv | awk 'NR > 2 {print NF; exit}')
        last_c=\$((\$numb-1))
        gc_c=\$((\$numb+1))

        #get GC content of called nucleosomes
        awk -v OFS='\t' 'NR > 2 {for(i=2;i<=NF;++i)printf \$i""FS; print}' subNucs_readCounts.csv | cut -f 2- \
        | bedtools nuc -fi $params.genome -bed - | cut -f 1-\$last_c,\$gc_c \
        | awk -v OFS='\t' 'NR > 1 {print}'  > pre_subNucs_readCounts_wGC.csv

        #add nucID
        awk -v OFS='\t' 'NR > 2 {print \$1}' subNucs_readCounts.csv \
        | paste pre_subNucs_readCounts_wGC.csv - >pre_subNucs_readCounts_wGC_ext.csv

        #prepare header
        awk -v OFS='\t' 'FNR == 2 {print}' subNucs_readCounts.csv | cut -f 2-\$numb > header_featureCounts.csv
        echo "GC_cont\tnucID" | paste header_featureCounts.csv - > header_all.csv

        #add header
        cat header_all.csv  pre_subNucs_readCounts_wGC_ext.csv > subNucs_readCounts_wGC.csv

        """
}
