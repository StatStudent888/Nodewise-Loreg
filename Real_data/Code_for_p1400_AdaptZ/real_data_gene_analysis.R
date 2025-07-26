load(file="RData/L0_venn_plot.RData")
load(file="RData/L0_method3_venn_plot.RData")
load(file="RData/L1_venn_plot.RData")
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


library(stringr)
x1 <- setdiff(L1,L0)
probe_L1 <- str_split_fixed(x1,"-",n=2)
probe_L1_partial <- {}
for(i in 1:nrow(probe_L1)){
  probe_L1_partial_k <- NS_partial_correlation[which(rownames(NS_partial_correlation)==probe_L1[i,1]),which(colnames(NS_partial_correlation)==probe_L1[i,2])]
  probe_L1_partial <- c(probe_L1_partial,probe_L1_partial_k)
}
gene_L1 <- data.frame(matrix(NA,nrow=length(x1),ncol=0))
gene_L1$probe <- x1
gene_L1$partial <- probe_L1_partial
gene_L1

x2 <- setdiff(CLIME,L0)
probe_CLIME <- str_split_fixed(x2,"-",n=2)
probe_CLIME_partial <- {}
for(i in 1:nrow(probe_CLIME)){
  probe_CLIME_partial_k <- CLIME_partial_correlation[which(rownames(CLIME_partial_correlation)==probe_CLIME[i,1]),which(colnames(CLIME_partial_correlation)==probe_CLIME[i,2])]
  probe_CLIME_partial <- c(probe_CLIME_partial,probe_CLIME_partial_k)
}
gene_CLIME <- data.frame(matrix(NA,nrow=length(x2),ncol=0))
gene_CLIME$probe <- x2
gene_CLIME$partial <- probe_CLIME_partial
gene_CLIME

x3 <- setdiff(Glasso,L0)
probe_Glasso <- str_split_fixed(x3,"-",n=2)
probe_Glasso_partial <- {}
for(i in 1:nrow(probe_Glasso)){
  probe_Glasso_partial_k <- Glasso_partial_correlation[which(rownames(Glasso_partial_correlation)==probe_Glasso[i,1]),which(colnames(Glasso_partial_correlation)==probe_Glasso[i,2])]
  probe_Glasso_partial <- c(probe_Glasso_partial,probe_Glasso_partial_k)
}
gene_Glasso <- data.frame(matrix(NA,nrow=length(x3),ncol=0))
gene_Glasso$probe <- x3
gene_Glasso$partial <- probe_Glasso_partial
gene_Glasso

probe_total <- unique(c(x1,x2,x3))
gene_total <- data.frame(matrix(NA,nrow=length(probe_total),ncol=0))
gene_total$probe <- probe_total
for(i in 1:nrow(gene_total)){
  index <- which(gene_L1$probe==gene_total$probe[i])
  if(length(index)==0){
    gene_total$L1[i] <- 0
  }
  else{gene_total$L1[i] <- gene_L1$partial[index]}
}

for(i in 1:nrow(gene_total)){
  index <- which(gene_CLIME$probe==gene_total$probe[i])
  if(length(index)==0){
    gene_total$CLIME[i] <- 0
  }
  else{gene_total$CLIME[i] <- gene_CLIME$partial[index]}
}

for(i in 1:nrow(gene_total)){
  index <- which(gene_Glasso$probe==gene_total$probe[i])
  if(length(index)==0){
    gene_total$Glasso[i] <- 0
  }
  else{gene_total$Glasso[i] <- gene_Glasso$partial[index]}
}

for(i in 1:nrow(gene_total)){
  gene_total$max[i] <- max(gene_total$L1[i],gene_total$CLIME[i],gene_total$Glasso[i])
}

for(i in 1:nrow(gene_total)){
  gene_total$median[i] <- median(gene_total$L1[i],gene_total$CLIME[i],gene_total$Glasso[i])
}

length(which(gene_total$max<=0.1))/nrow(gene_total)
length(which(gene_total$median<=0.1))/nrow(gene_total)

length(which(gene_total$max>0.1))
length(which(gene_total$median>0.1)) #54


L0_omega4_partial_correlation <- sdar_T_star1
for(i in 1:p){
  for(j in 1:p){
    if(i!=j){
      L0_omega4_partial_correlation[i,j] <- abs(-sdar_T_star1[i,j]/(sqrt(sdar_T_star1[i,i]*sdar_T_star1[j,j])))
    }
    else{L0_omega4_partial_correlation[i,j] <- 0}
  }
}


gene_54 <- gene_total[which(gene_total$median>0.1),]
gene_54
gene_54_probesplit <- str_split_fixed(gene_54$probe,"-",n=2) 
for(i in 1:nrow(gene_54)){
  gene_54$L0_omega4_partial[i] <- L0_omega4_partial_correlation[which(rownames(L0_omega4_partial_correlation)==gene_54_probesplit[i,1]),which(colnames(L0_omega4_partial_correlation)==gene_54_probesplit[i,2])]
}
length(which(gene_54$L0_omega4_partial!=0)) #3


