#Genetic_data_processing document includes original data, R code Genetic_data_processing for data pre-processing and screening for differentially expressed genes and Log2 document contains data after processing.

#Code_for_p300_BH document includes 5 R files and 1 RData document.
  #l0_desparsified_real_data_lfdr.R, l0_method3_real_data_lfdr.R, l0_unbias_real_data_lfdr.R compute all 5 estimators of Nodewise Loreg algorithm in real   
    data with 300 selected genes based on BH and return their classification performance and sparsity.
  #l0_unbias_real_data_network-graph.R draw network-graph based on estimated absolute partial correlations and estimated absolute partial 
    correlations thresholded by 0.1.
  #Gene_analysis_and_VennDiagram.R draw Venndiagram based on estimated absolute partial correlations and estimated absolute partial correlations 
    thresholded by 0.1. Moreover, this code gives all results of the analysis of direct connectivity in our paper.
  #RData document includes all Rdata used in Gene_analysis_and_VennDiagram.R.

#Code_for_p1400_AdaptZ document includes 8 R files and 1 RData document.
  #l0_desparsified_real_data_p1400.R, l0_method3_real_data_p1400.R, l0_unbias_real_data_p1400.R compute all 5 estimators of Nodewise Loreg  
    algorithm in real data with 1400 selected genes based on AdaptZ and return their classification performance and sparsity.
  #l0_unbias_real_data_network-graph_p1400.R draw network-graph based on estimated absolute partial correlations and estimated absolute partial 
    correlations thresholded by 0.1.
  #l0_real_data_p1400_time.R computes running time of Nodewise Loreg estimator under one replication and saves the result and l0_real_data_p1400 
    _time_combine.R combines the results of running time under 100 replication.
  #real_data_gene_analysis.R and VennDiagram.R draw Venndiagram based on estimated absolute partial correlations and estimated absolute partial  
    correlations thresholded by 0.1 and give all results of the analysis of direct connectivity in our paper.
  #RData document includes all Rdata used in real_data_gene_analysis.R and VennDiagram.R


