/*
* Show settings at the beginning
*/

def settings() {
        println ''
         log.info """\

                  nucMACC   P I P E L I N E
                  =============================
                  Path variables used in analysis
                  csvInput :    ${params.csvInput}
                  outDir   :    ${params.outDir}
                  genomeIdx:    ${params.genomeIdx}

                  General options


                  """.stripIndent()

        println ''
}
