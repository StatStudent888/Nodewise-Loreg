#The document FDR-R-Code includes code of AdaptZ (Sun and Cai, 2007) to control FDR.
#The document l0bnb-master includes python code for running L0BnB algorithm (Hazimeh et al., 2022).
#The document subset_selection_with_shrinkage-master includes python code for running MIO algorithm (Mazumder et al., 2023).
#abnb-end.R includes a function that runs Nodewise L0BnB algorithm and returns the optimal estimation of coefficient using HBIC criterion.
#amio-end.R includes a function that runs Nodewise MIO algorithm and returns the optimal estimation of coefficient using HBIC criterion.
#ans-end.R includes a function that runs Nodewise Lasso (L1) algorithm and returns the optimal estimation of coefficient using HBIC criterion.
#asdar.cpp includes functions that runs Nodewise Loreg (L0) algorithm and returns the optimal estimation of coefficient using HBIC criterion. There is no need to input the normalized version of design matrix X because the code will automatically normalize design matrix X.
#asdar_youhua_Eigen.cpp includes the corresponding functions of asdar written in the Eigen version, which is faster than armadillo, and is used to compute the running time for Nodewise Loreg. 
#compute_desparsified.cpp includes functions to compute the true length of confidence interval for Nodewise Loreg under both Gaussian and sub-Gaussian design. 
#evaluation.R includes a function that returns the support recovery performance of estimators on off-diagonal entries.
#matrix_generator.R includes functions that generate symmetric matrix or asymmetric matrix after inputting a vector that contains diagonal elements of a matrix and a (p-1)*p matrix that does not contain diagonal elements.
#omega_generator.R includes a function that generate Band graph model.
#real_evaluation.R includes a function that returns classification performance used in real data.
#score.R includes a function that returns LDA score used in real data.
#useful_function.R includes two functions which normalize each column of a matrix to have unit L2 norm the matrix and do matrix multiplication.