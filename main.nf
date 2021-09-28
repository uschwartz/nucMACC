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

// load modules
//fastqc
include{fastqc; multiqc} from './modules/raw_qc'
//alignment to ref genome
include{alignment} from './modules/align'
//qualimap after alignment
include{qualimap} from './modules/qualimap'
// filtering sizes using alignmentSieve
include{sieve_mono; sieve_sub} from './modules/alignmentsieve'
// prepare for DANPOS
include{pool} from './modules/prepareDANPOS'
// DANPOS run
include{danpos_mono} from './modules/DANPOS'


workflow{
  fastqc(sampleSingle_ch)
  alignment(samplePair_ch)
  qualimap(alignment.out[1])
  sieve_mono(alignment.out[1])
  sieve_sub(alignment.out[1])
  multiqc(fastqc.out[0].mix(alignment.out[0]).mix(qualimap.out).collect())
  pool(sieve_mono.out[1].map{name,bam -> file(bam)}.collect())
  danpos_mono(pool.out[0], pool.out[1])
}
