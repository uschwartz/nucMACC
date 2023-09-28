//load modules
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
//prepare for DANPOS
include{extract_chrSizes} from '../modules/prepareDANPOS'
//DANPOS run
include{danpos_mono; danpos_sub} from '../modules/DANPOS'
// get read count per nucleosome
include{featureCounts_mono; featureCounts_sub} from '../modules/featureCounts'
//deeptools TSS
include{TSS_profile_mono;TSS_profile_plot_mono;TSS_profile_sub;TSS_profile_plot_sub} from '../modules/get_TSS_profile'
//TSS Profile monoNucs
include{make_TSS_plots_monoNucs;make_TSS_plots_subNucs } from '../modules/make_TSS_plots'
//bigwigSummary for correlationHeatmap and PCA
include{bigwigSummary as bigwigSummary_all; bigwigSummary as bigwigSummary_monoNucs; bigwigSummary as bigwigSummary_subNucs} from '../modules/bigwigSummary'
//correlation heatmap
include{plot_correlation_Heatmap as plot_correlation_Heatmap_all; plot_correlation_Heatmap as plot_correlation_Heatmap_monoNucs; plot_correlation_Heatmap as plot_correlation_Heatmap_subNucs} from '../modules/correlationHeatmap'
// PCA
include{plot_PCA as plot_PCA_all; plot_PCA as plot_PCA_monoNucs; plot_PCA as plot_PCA_subNucs} from '../modules/PCA'


workflow MNaseQC{

  take:
  sampleSingle_ch
  samplePair_ch
  samples_conc

  main:
  fastqc(sampleSingle_ch)
  alignment(samplePair_ch)
  qualimap(alignment.out[1])


  // monoNucs
  sieve_mono(alignment.out[1])
  extract_chrSizes(alignment.out[1].collect())
  danpos_mono(sieve_mono.out[1], extract_chrSizes.out[0])

  //subNucs
  sieve_sub(alignment.out[1])
  danpos_sub(sieve_sub.out[1], extract_chrSizes.out[0])

  //QualityCheck
  multiqc(fastqc.out[0].mix(alignment.out[0]).mix(qualimap.out[0]).collect())
  InsertSize_Histogram(qualimap.out[0].collect())

  //FragmentStatistics
  statistics_read(sieve_mono.out[0].join(sieve_sub.out[0]).join(fastqc.out[2]).join(alignment.out[2]).join(qualimap.out[1]))
  statistics_plot(statistics_read.out[0].collect())

  //TSS_Profile_mono
  if(params.TSS){
  TSS_profile_mono(danpos_mono.out[0].collect())
  TSS_profile_plot_mono(TSS_profile_mono.out)
  make_TSS_plots_monoNucs(TSS_profile_plot_mono.out)

  //TSS_Profile_sub
  TSS_profile_sub(danpos_sub.out[0].collect())
  TSS_profile_plot_sub(TSS_profile_sub.out)
  make_TSS_plots_subNucs(TSS_profile_plot_sub.out) }

  //bigwigSummary
  bigwigSummary_all(danpos_mono.out[0].mix(danpos_sub.out[0]).collect(), Channel.value('all'))
  bigwigSummary_monoNucs(danpos_mono.out[0].collect(), Channel.value('monoNucs'))
  bigwigSummary_subNucs(danpos_sub.out[0].collect(), Channel.value('subNucs'))

  //correlation heatmap with Spearman correlation
  plot_correlation_Heatmap_all(bigwigSummary_all.out[0], Channel.value('all'))
  plot_correlation_Heatmap_monoNucs(bigwigSummary_monoNucs.out[0], Channel.value('monoNucs'))
  plot_correlation_Heatmap_subNucs(bigwigSummary_subNucs.out[0], Channel.value('subNucs'))

  //plot PCA
  plot_PCA_all(bigwigSummary_all.out[0], Channel.value('all'))
  plot_PCA_monoNucs(bigwigSummary_monoNucs.out[0], Channel.value('monoNucs'))
  plot_PCA_subNucs(bigwigSummary_subNucs.out[0], Channel.value('subNucs'))
}
