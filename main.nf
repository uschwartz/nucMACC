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
// forward reads
Channel
    .fromPath(params.csvInput)
    .splitCsv(header:true)
    .map{ row -> tuple(row.Sample_Name,
      file(params.project.concat(row.path_fwdReads)))}
    .set{samples_fwd_ch}

// reverse reads
Channel
      .fromPath(params.csvInput)
      .splitCsv(header:true)
      .map{ row -> tuple(row.Sample_Name,
         file(params.project.concat(row.path_revReads)))}
      .set{samples_rev_ch}

//Channel for fastqc
samples_fwd_ch.mix(samples_rev_ch).set{sampleSingle_ch}
//Channel for alignment
samples_fwd_ch.join(samples_rev_ch).set{samplePair_ch}

//read MNase concentration
Channel
      .fromPath(params.csvInput)
      .splitCsv(header:true)
      .map{ row -> tuple(row.MNase_U,row.Sample_Name)}
      .set{samples_conc}

//get lowest MNase digest
samples_conc.map{conc,sample -> conc}.min().set{min_conc}

// get sample with lowest MNase digest
min_conc.join(samples_conc)
.map{conc,sample -> sample}
.set{min_conc_sample}


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
include{danpos_mono; danpos_sub} from './modules/DANPOS'
// convert DANPOS output
include{convert2saf_mono; convert2saf_sub} from './modules/convert2saf'
// get read count per nucleosome
include{featureCounts_mono; featureCounts_sub} from './modules/featureCounts'


workflow{
  fastqc(sampleSingle_ch)
  alignment(samplePair_ch)
  qualimap(alignment.out[1])

  // monoNucs
  sieve_mono(alignment.out[1])
  pool(sieve_mono.out[1].map{name,bam -> file(bam)}.collect())
  danpos_mono(sieve_mono.out[1].mix(pool.out[0]), pool.out[1])
  convert2saf_mono(danpos_mono.out[1].join(pool.out[0]))
  featureCounts_mono(convert2saf_mono.out[1], sieve_mono.out[1].map{name,bam -> file(bam)}.collect())

  //subNucs
  sieve_sub(alignment.out[1])
  danpos_sub(sieve_sub.out[1], pool.out[1])
  convert2saf_sub(danpos_sub.out[1].join(min_conc_sample).join(sieve_sub.out[1]))
  featureCounts_sub(convert2saf_sub.out[1], sieve_sub.out[1].map{name,bam -> file(bam)}.collect())

  multiqc(fastqc.out[0].mix(alignment.out[0]).mix(qualimap.out).collect())
}
