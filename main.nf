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
if(params.test){
  if (params.bamEntry){
          // sample mono
          Channel
              .fromPath(params.csvInput)
              .splitCsv(header:true)
              .map{ row -> tuple(row.Sample_Name,file(row.path_mono))}
              .set{bamEntry_mono}
          // sample sub
          Channel
              .fromPath(params.csvInput)
              .splitCsv(header:true)
              .map{ row -> tuple(row.Sample_Name,file(row.path_sub))}
              .set{bamEntry_sub}
          }
  else {
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
      }
}

else if (params.bamEntry == true){
        // sample mono
        Channel
            .fromPath(params.csvInput)
            .splitCsv(header:true)
            .map{ row -> tuple(row.Sample_Name,file(row.path_mono))}
            .set{bamEntry_mono}
        // sample sub
        Channel
            .fromPath(params.csvInput)
            .splitCsv(header:true)
            .map{ row -> tuple(row.Sample_Name,file(row.path_sub))}
            .set{bamEntry_sub}
        println "BamEntry csv part"  }


else {
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
}

if(params.bamEntry==false) {
  //Channel for fastqc
  samples_fwd_ch.mix(samples_rev_ch).set{sampleSingle_ch}
  //Channel for alignment
  samples_fwd_ch.join(samples_rev_ch).set{samplePair_ch}
}


//read MNase concentration
Channel
      .fromPath(params.csvInput)
      .splitCsv(header:true)
      .map{ row -> tuple(row.MNase_U.toDouble(),row.Sample_Name)}
      .set{samples_conc}


// load workflows
// generate profiles
include{MNaseQC} from './workflows/MNaseQC'
include{sub_bamEntry; sub_FASTQ_entry; common_nucMACC} from './workflows/nucMACC'


// Check mandatory parameters
if (params.csvInput) { ch_csv = file(params.csvInput) } else { exit 1, 'Input samplesheet not found!' }

if (params.TSS) {ch_TSS = file (params.TSS)}
      if(params.TSS){
        if (ch_TSS.isEmpty()) { exit 1, 'TSS file not found!'}
        }

if (params.blacklist) {ch_blacklist = file (params.blacklist)}
      if(params.blacklist){
        if (ch_blacklist.isEmpty()) { exit 1, 'Blacklist file not found!'}
        }
  if(params.analysis =='MNaseQC'){
    if (params.genomeIdx) { ch_idx = file(params.genomeIdx).parent } 
    if (ch_idx.isEmpty()) { exit 1, 'Folder containing bowtie2 indices not found!'}
  }

  if(params.analysis=='nucMACC'){
    if(params.bamEntry == true){

    if (params.genome) { ch_genome = file(params.genome)}
    if (ch_genome.isEmpty()) { exit 1, 'Genome fasta not found!'}

    }
    else {
      if (params.genomeIdx) { ch_idx = file(params.genomeIdx).parent } 
    if (ch_idx.isEmpty()) { exit 1, 'Folder containing bowtie2 indices not found!'}
      if (params.genome) { ch_genome = file(params.genome)}
    if (ch_genome.isEmpty()) { exit 1, 'Genome fasta not found!'}}
    }
  
  

workflow{
        if(params.analysis=='MNaseQC'){
                MNaseQC(sampleSingle_ch,samplePair_ch,samples_conc)
        }
        if(params.analysis=='nucMACC'){
    


                if (params.bamEntry == true) {
                  sub_bamEntry(bamEntry_mono,bamEntry_sub)
                  common_nucMACC(sub_bamEntry.out[0], sub_bamEntry.out[1], samples_conc)
                  }
                else {
                  sub_FASTQ_entry(sampleSingle_ch,samplePair_ch,samples_conc)
                  common_nucMACC(sub_FASTQ_entry.out[0], sub_FASTQ_entry.out[1], samples_conc)
                }
        }
}
