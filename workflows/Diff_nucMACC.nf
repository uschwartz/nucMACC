//load Modules
// reference Map
include{reference_map_mono; reference_map_sub} from '../modules/reference_map'
//get counts
include{featureCounts_diff_mono; featureCounts_diff_sub} from '../modules/featureCounts_diff'
//include Diff analysis
include{diff_nucMACC_mono; diff_nucMACC_sub} from '../modules/diff_nucMACC'
//include score comparison
include{compare_nucMACC_mono; compare_nucMACC_sub} from '../modules/compare_nucMACC'


workflow Diff_nucMACC{
    take:
    bamEntry_mono
    bamEntry_sub
   
    main:
    //monoNucs
    //get Reference map of extreme Nucs
    reference_map_mono(Channel.fromPath(params.csvInput))
    //get read count
    featureCounts_diff_mono(reference_map_mono.out[1], bamEntry_mono.map{Sample_Name,path_mono -> file(path_mono)}.collect(),bamEntry_mono.map{Sample_Name,path_mono -> Sample_Name}.collect())
    //diff Analysis
    diff_nucMACC_mono(Channel.fromPath(params.csvInput),featureCounts_diff_mono.out[0],featureCounts_diff_mono.out[1])
    //compare Scores
    compare_nucMACC_mono(Channel.fromPath(params.csvInput), diff_nucMACC_mono.out[2], diff_nucMACC_mono.out[1])



    //subNucs
    //get Reference map of extreme Nucs
    reference_map_sub(Channel.fromPath(params.csvInput))
    //get read count
    featureCounts_diff_sub(reference_map_sub.out[1], bamEntry_sub.map{Sample_Name,path_sub -> file(path_sub)}.collect(),bamEntry_sub.map{Sample_Name,path_sub -> Sample_Name}.collect())
    //diff Analysis
    diff_nucMACC_sub(Channel.fromPath(params.csvInput),featureCounts_diff_sub.out[0],featureCounts_diff_sub.out[1])
    //compare scores
    compare_nucMACC_sub(Channel.fromPath(params.csvInput), diff_nucMACC_sub.out[2], diff_nucMACC_sub.out[1])
}
