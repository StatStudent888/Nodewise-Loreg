#A function that normalizes each column of a matrix to have unit L2 norm
cppFunction('
  List normalizeMatrix(NumericMatrix X) {
    int n = X.nrow(); 
    int m = X.ncol(); 
    NumericMatrix normalizedX(n, m); 
    NumericVector norms(m); 
    
    for (int j = 0; j < m; j++) {
      double colNorm = 0.0;
      for (int i = 0; i < n; i++) {
        colNorm += X(i, j) * X(i, j);
      }
      norms[j] = 1 / sqrt(colNorm);
    }
    
    for (int j = 0; j < m; j++) {
      for (int i = 0; i < n; i++) {
        normalizedX(i, j) = X(i, j) * norms[j];
      }
    }
    
    return List::create(
      Named("normalizedX") = normalizedX,  
      Named("norms") = norms              
    );
  }
')

#A function for fast matrix multiplication
cppFunction('
  NumericVector dotMultiply(NumericVector dx, NumericVector a) {
    int n = dx.size();  
    NumericVector result(n);  
    
    for (int i = 0; i < n; i++) {
      result[i] = dx[i] * a[i];  
    }
    
    return result;
  }
')




