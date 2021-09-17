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
         include{settings} from './modules/setting'
         settings()
 }

 //                       help message
 // Show help message
 if (params.help) {
     include{helpMessage} from './modules/help'
     helpMessage()
     exit 0
 }

 //                      workflow



Channel
    .fromPath(params.csvInput)
    .splitCsv(header:true)
    .map{ row -> tuple(row.Sample_Name,
      file(params.project.concat(row.path_fwdReads)))}
    .set{samples_fwd_ch}


Channel
      .fromPath(params.csvInput)
      .splitCsv(header:true)
      .map{ row -> tuple(row.Sample_Name,
         file(params.project.concat(row.path_revReads)))}
      .set{samples_rev_ch}

samples_fwd_ch.join(samples_rev_ch).set{samplePair_ch}
samples_fwd_ch.mix(samples_rev_ch).set{sampleSingle_ch}


process fastqc{
  echo true

  input:
  tuple val(sampleID), file(read)

  output:
  stdout

  script:
  """
  fastqc $read
  """
}

workflow{
  fastqc(sampleSingle_ch)
}
