process featureCounts_mono{
        label 'big'

        publishDir "${params.outDir}/RUN/04_NUCS_READ_COUNTS", mode: 'copy'


        input:
        file(saf)
        file(monoNucs)

        output:
        file("monoNucs_readCounts_wGC.csv")
        file("monoNucs_readCounts_all.csv.summary")

        script:
        """
        featureCounts  -F SAF -a $saf \
        -o monoNucs_readCounts_all.csv \
        --fracOverlap 0.7 \
        -T $task.cpus -p -B -O --largestOverlap \
        $monoNucs

        

        #filter by readcounts
        awk -F'\t' 'NR > 2 { sum = 0; for (i = 7; i <= NF; i++) {sum += \$i;} if (sum >= $params.rawfilt_mono) { print; } }' monoNucs_readCounts_all.csv > monoNucs_readCounts_pre.csv
        
        #store first 2 lines as they get lost
        awk 'NR <= 2' monoNucs_readCounts_all.csv > monoNucs_readCounts_head.csv

        #readd first two lines
        cat monoNucs_readCounts_head.csv monoNucs_readCounts_pre.csv > monoNucs_readCounts.csv
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

        publishDir "${params.outDir}/RUN/04_NUCS_READ_COUNTS", mode: 'copy'


        input:
        file(saf)
        file(subNucs)
        val(lowestcond)

        output:
        file("subNucs_readCounts_wGC.csv")
        file("subNucs_readCounts_all.csv.summary")

        script:
        """
        featureCounts  -F SAF -a $saf \
        -o subNucs_readCounts_all.csv \
        --fracOverlap 0.7 \
        -T $task.cpus -p -B -O --largestOverlap \
        $subNucs


        export lowest_cond_col=\$(awk -F'\t' 'NR == 2 { for (i = 1; i <= NF; i++) { if (\$i ~ /$lowestcond/) { j = i; break; } } } END { print j; }' subNucs_readCounts_all.csv)
        awk -v col="\$lowest_cond_col" -v rawfilt="$params.rawfilt_sub" -F'\t' 'NR <= 2 || (\$col >= rawfilt) {print;}' subNucs_readCounts_all.csv > subNucs_readCounts.csv


        
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
        awk -v OFS='\t' 'FNR == 2 {print}' subNucs_readCounts_all.csv | cut -f 2-\$numb > header_featureCounts.csv
        echo "GC_cont\tnucID" | paste header_featureCounts.csv - > header_all.csv

        #add header
        cat header_all.csv  pre_subNucs_readCounts_wGC_ext.csv > subNucs_readCounts_wGC.csv

        """
}
