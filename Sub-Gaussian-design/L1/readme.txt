#One-replication document includes 3 R code:
    #L1_desparsified_subGaussian_normality.R computes AvgLength, CovRate, AbsAvgZ, SDZ of desparsified Nodewise Lasso estimator under  
      sub-Gaussian design under one replication and saves the result.
    #L1_desparsified_subGaussian_normalplot.R computes Z-score of desparsified Nodewise Lasso estimator under sub-Gaussian design under  
      one replication and saves the result.
    #L1_desparsified_subGaussian_one-two-stage_lfdr.R computes all 5 estimators of Nodewise Lasso algorithm and returns their support                      
      recovery performance and matrix norm losses under sub-Gaussian design under one replication and saves the result. 

#Combine document includes 3 corresponding R code that combine the results under 400 replications for normal plot and 100 replications for    
   others based on the preserved Rdata using code in One-replication document.

