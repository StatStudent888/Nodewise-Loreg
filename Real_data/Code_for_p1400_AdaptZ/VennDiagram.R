load(file="RData/L0_venn_plot.RData")
load(file="RData/L1_venn_plot.RData")
load(file="RData/L1_method3_venn_plot.RData")
load(file="RData/CLIME_venn_plot.RData")
load(file="RData/Glasso_venn_plot.RData")

p <- 1400

L0 <- {}
L1 <- {}
CLIME <- {}
Glasso <- {}
for(i in 2:p){
  for(j in 1:(i-1)){
    if(sdar_unbias_partial_correlation_star[i,j]!=0){
      L0_k <- paste(rownames(sdar_unbias_partial_correlation_star)[i],colnames(sdar_unbias_partial_correlation_star)[j],sep="-")
      L0 <- c(L0,L0_k)
    }
  }
}

for(i in 2:p){
  for(j in 1:(i-1)){
    if(NS_partial_correlation[i,j]!=0){
      L1_k <- paste(rownames(NS_partial_correlation)[i],colnames(NS_partial_correlation)[j],sep="-")
      L1 <- c(L1,L1_k)
    }
  }
}

for(i in 2:p){
  for(j in 1:(i-1)){
    if(CLIME_partial_correlation[i,j]!=0){
      CLIME_k <- paste(rownames(CLIME_partial_correlation)[i],colnames(CLIME_partial_correlation)[j],sep="-")
      CLIME <- c(CLIME,CLIME_k)
    }
  }
}

for(i in 2:p){
  for(j in 1:(i-1)){
    if(Glasso_partial_correlation[i,j]!=0){
      Glasso_k <- paste(rownames(Glasso_partial_correlation)[i],colnames(Glasso_partial_correlation)[j],sep="-")
      Glasso <- c(Glasso,Glasso_k)
    }
  }
}

x1 <- setdiff(L0,L1)
x2 <- setdiff(x1,CLIME)
x3 <- setdiff(x2,Glasso)
x3
library(stringr)
probe_id <- str_split_fixed(x3,"-",n=2)
gene <- data.frame(matrix(NA,nrow=length(x3),ncol=3))
colnames(gene) <- c("probe","partial","L1_Omega4_partial")
gene$probe <- probe_id
gene

for(i in 1:length(x3)){
  gene$partial[i] <- sdar_unbias_partial_correlation_star[which(rownames(sdar_unbias_partial_correlation_star)==probe_id[i,1]),which(colnames(sdar_unbias_partial_correlation_star)==probe_id[i,2])]
}
gene

NS_method3_partial_correlation <- NS_T_star1
for(i in 1:p){
  for(j in 1:p){
    if(i!=j){
      NS_method3_partial_correlation[i,j] <- -NS_T_star1[i,j]/(sqrt(NS_T_star1[i,i]*NS_T_star1[j,j]))
    }
    else{NS_method3_partial_correlation[i,j] <- 0}
  }
}
for(i in 1:length(x3)){
  gene$L1_Omega4_partial[i] <- NS_method3_partial_correlation[which(rownames(NS_method3_partial_correlation)==probe_id[i,1]),which(colnames(NS_method3_partial_correlation)==probe_id[i,2])]
}

gene
nrow(gene[which(gene$L1_Omega4_partial!=0),]) #310
gene[which(gene$L1_Omega4_partial==0),]
nrow(gene[which(gene$L1_Omega4_partial==0),]) #82
gene_82 <- cbind(probe_id[which(gene$L1_Omega4_partial==0),],gene[which(gene$L1_Omega4_partial==0),])
gene_82_order <- gene_82[order(abs(gene_82$partial),decreasing=TRUE),]
gene_82_order
length(which(gene$L1_Omega4_partial==0))


load(file="RData/L0_venn_plot.RData")
load(file="RData/L1_venn_plot.RData")
load(file="RData/CLIME_venn_plot.RData")
load(file="RData/Glasso_venn_plot.RData")

p <- 1400

L0_1 <- {}
L1_1 <- {}
CLIME_1 <- {}
Glasso_1 <- {}
for(i in 2:p){
  for(j in 1:(i-1)){
    if(sdar_unbias_partial_correlation_star[i,j]!=0){
      L0_1k <- paste(rownames(sdar_unbias_partial_correlation_star)[i],colnames(sdar_unbias_partial_correlation_star)[j],sep="-")
      L0_1 <- c(L0_1,L0_1k)
    }
  }
}

for(i in 2:p){
  for(j in 1:(i-1)){
    if(NS_partial_correlation[i,j]!=0){
      L1_1k <- paste(rownames(NS_partial_correlation)[i],colnames(NS_partial_correlation)[j],sep="-")
      L1_1 <- c(L1_1,L1_1k)
    }
  }
}

for(i in 2:p){
  for(j in 1:(i-1)){
    if(CLIME_partial_correlation[i,j]!=0){
      CLIME_1k <- paste(rownames(CLIME_partial_correlation)[i],colnames(CLIME_partial_correlation)[j],sep="-")
      CLIME_1 <- c(CLIME_1,CLIME_1k)
    }
  }
}

for(i in 2:p){
  for(j in 1:(i-1)){
    if(Glasso_partial_correlation[i,j]!=0){
      Glasso_1k <- paste(rownames(Glasso_partial_correlation)[i],colnames(Glasso_partial_correlation)[j],sep="-")
      Glasso_1 <- c(Glasso_1,Glasso_1k)
    }
  }
}

x <- list(
  A = L0_1,
  B = L1_1,
  C = CLIME_1,
  D = Glasso_1
)

display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

display_venn(
  x,
  category.names = c(expression(L[0]) , expression(L[1]) , "CLIME", "GLasso"),
  fill = c("#999999", "#E69F00", "#56B4E9", "#009E73")
)




