#!/usr/bin/env nextflow

 /*
 ===============================================================================
                      nextflow based nucMACC pipeline
 ===============================================================================
Authors:
Uwe Schwartz <uwe.schwartz@ur.de>
 -------------------------------------------------------------------------------
 */

nextflow.enable.dsl = 2

 //                           show settings
 if (!params.help) {
 println ''
  log.info """\

           nucMACC   P I P E L I N E
           =============================
           Path variables used in analysis
           csvInput : ${params.csvInput}


           General options


           """.stripIndent()

 println ''
 }


 //                       help message

 /*
 * Help message
 */
 def helpMessage() {
     println ''
     log.info """
     nucMACC   P I P E L I N E
     =============================
     Usage:

     nextflow run uschwartz/nucMACC --csvInput 'path2csvFile'

     Mandatory arguments:
       --csvInput           [string] Path to comma-separated file containing information about the samples in the experiment (see ./toyData/input.csv as example) template provided in ./input_template.csv

      """.stripIndent()
      println ''
 }

 // Show help message
 if (params.help) {
     helpMessage()
     exit 0
 }
