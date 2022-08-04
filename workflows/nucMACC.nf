// load modules
//fastqc
include{fastqc; multiqc} from '../modules/raw_qc'
//alignment to ref genome
include{alignment} from '../modules/align'
//qualimap after alignment
include{qualimap} from '../modules/qualimap'
//InsertSize_Histogram
include{InsertSize_Histogram} from '../modules/InsertSize_Histogram'
//FragmentStatistics
include{statistics_read; statistics_plot} from '../modules/fragment_statistics'
// filtering sizes using alignmentSieve
include{sieve_mono; sieve_sub} from '../modules/alignmentsieve'
// prepare for DANPOS
include{pool} from '../modules/prepareDANPOS'
// DANPOS run
include{danpos_mono; danpos_sub} from '../modules/DANPOS'
// convert DANPOS output
include{convert2saf_mono; convert2saf_sub} from '../modules/convert2saf'
// get read count per nucleosome
include{featureCounts_mono; featureCounts_sub} from '../modules/featureCounts'
// get nucMACC_scores
include{nucMACC_scores;sub_nucMACC_scores} from '../modules/nucMACC_scores'
//deeptools TSS
include{TSS_profile_mono;TSS_profile_plot_mono;TSS_profile_sub;TSS_profile_plot_sub} from '../modules/get_TSS_profile'
//TSS Profile monoNucs
include{make_TSS_plots_monoNucs;make_TSS_plots_subNucs } from '../modules/make_TSS_plots'

// bamEntry module
include{mergeBam_mono; mergeBam_sub} from '../modules/bamEntry'

workflow sub_bamEntry{

    take:
    bamEntry_mono
    bamEntry_sub

    main:
    mergeBam_mono(bamEntry_mono.groupTuple())
    mergeBam_sub(bamEntry_sub.groupTuple())

    emit:
    mergeBam_mono.out
    mergeBam_sub.out

  }

  workflow sub_FASTQ_entry{

    take:
    sampleSingle_ch
    samplePair_ch
    samples_conc

    main:
    fastqc(sampleSingle_ch)
    alignment(samplePair_ch)
    qualimap(alignment.out[1])

    //monoNuc
    sieve_mono(alignment.out[1])

    //subNucs
    sieve_sub(alignment.out[1])

    //QualityCheck
    multiqc(fastqc.out[0].mix(alignment.out[0]).mix(qualimap.out[0]).collect())
    InsertSize_Histogram(qualimap.out[0].collect())

    //FragmentStatistics
    statistics_read(sieve_mono.out[0].join(sieve_sub.out[0]).join(fastqc.out[2]).join(alignment.out[2]).join(qualimap.out[1]))
    statistics_plot(statistics_read.out[0].collect())

    emit:
    sieve_mono.out[1]
    sieve_sub.out[1]

  }

workflow common_nucMACC{

  take:
  sieve_mono
  sieve_sub
  samples_conc

  main:
  //get lowest MNase digest
  samples_conc.map{conc,sample -> conc}.min().set{min_conc}

  // get sample with lowest MNase digest
  min_conc.join(samples_conc)
  .map{conc,sample -> sample}
  .set{min_conc_sample}

  // monoNucs
  pool(sieve_mono.map{name,bam -> file(bam)}.collect())
  danpos_mono(sieve_mono.mix(pool.out[0]), pool.out[1])
  convert2saf_mono(danpos_mono.out[1].join(pool.out[0]))
  featureCounts_mono(convert2saf_mono.out[1], sieve_mono.map{name,bam -> file(bam)}.collect())

  // subNucs
  danpos_sub(sieve_sub, pool.out[1])
  convert2saf_sub(danpos_sub.out[1].join(min_conc_sample).join(sieve_sub))
  featureCounts_sub(convert2saf_sub.out[1], sieve_sub.map{name,bam -> file(bam)}.collect())

  //TSS_Profile_mono
  if(params.TSS){
  TSS_profile_mono(danpos_mono.out[0].collect())
  TSS_profile_plot_mono(TSS_profile_mono.out)
  make_TSS_plots_monoNucs(TSS_profile_plot_mono.out)
  //TSS_Profile_sub
  TSS_profile_sub(danpos_sub.out[0].collect())
  TSS_profile_plot_sub(TSS_profile_sub.out)
  make_TSS_plots_subNucs(TSS_profile_plot_sub.out)
  }

  //nucMACC scores
  nucMACC_scores(featureCounts_mono.out[0], Channel.fromPath(params.csvInput),featureCounts_mono.out[1])

  //subMACC scores
  sub_nucMACC_scores(featureCounts_sub.out[0], Channel.fromPath(params.csvInput),min_conc_sample, nucMACC_scores.out[2],featureCounts_mono.out, featureCounts_sub.out[1])

}
