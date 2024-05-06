#asdar-end.R includes a function that runs SDAR algorithm and returns the optimal estimation of coefficient using HBIC criterion. There is no need to input the normalized version of design matrix X because the code will automatically normalize design matrix X.
#evaluation.R includes a function that returns the support recovery performance of estimators on off-diagonal entries.
#matrix_generator.R includes functions that generate symmetric matrix or asymmetric matrix after inputting a vector that contains diagonal elements of a matrix and a (p-1)*p matrix that does not contain diagonal elements.
#omega_generator.R includes a function that generate Band graph model.
#real_evaluation.R includes a function that returns classification performance used in real data.
#score.R includes a function that returns LDA score used in real data.