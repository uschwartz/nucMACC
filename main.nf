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
// read csv file
if(params.analysis=='nucMACC'){
        // forward reads
        Channel
            .fromPath(params.csvInput)
            .splitCsv(header:true)
            .map{ row -> tuple(row.Sample_Name,file(row.path_fwdReads))}
            .set{samples_fwd_ch}
        // reverse reads
        Channel
              .fromPath(params.csvInput)
              .splitCsv(header:true)
              .map{ row -> tuple(row.Sample_Name,file(row.path_revReads))}
              .set{samples_rev_ch}
//Channel for fastqc
samples_fwd_ch.mix(samples_rev_ch).set{sampleSingle_ch}
//Channel for alignment
samples_fwd_ch.join(samples_rev_ch).set{samplePair_ch}

//read MNase concentration
Channel
      .fromPath(params.csvInput)
      .splitCsv(header:true)
      .map{ row -> tuple(row.MNase_U.toDouble(),row.Sample_Name)}
      .set{samples_conc}
}

// load workflows
// generate profiles
include{nucMACC} from './workflows/nucMACC'


workflow{
        if(params.analysis=='nucMACC'){
                profiler(sampleSingle_ch,samplePair_ch,samples_conc)
        }
        
}
