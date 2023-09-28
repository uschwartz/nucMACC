/*
* Show settings at the beginning
*/

def settings() {
        println ''
         log.info """\

                  nucMACC   P I P E L I N E
                  =============================
                  General options:
                  analysis         :  ${params.analysis}
                  test             :  ${params.test}
                  genomeSize       :  ${params.genomeSize}

                  Path variables used in analysis:
                  csvInput         :  ${params.csvInput}
                  outDir           :  ${params.outDir}
                  genomeIdx        :  ${params.genomeIdx}
                  genome           :  ${params.genome}

                  Additional options:
                  blacklist        :  ${params.blacklist}
                  TSS              :  ${params.TSS}
                  publishBam       :  ${params.publishBam}
                  publishBamFlt    :  ${params.publishBamFlt}

                  nucMACC specific options:
                  bamEntry         :  ${params.bamEntry}

                  MNaseQC specific options:
                  correlationMethod:  ${params.correlationMethod}
                  rawfilt_sub      :  ${params.rawfilt_sub}
                  rawfilt_mono     :  ${params.rawfilt_mono}


                  """.stripIndent()

        println ''
}
