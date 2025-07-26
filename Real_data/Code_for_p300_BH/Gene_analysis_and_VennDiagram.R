load(file="RData/L0_venn_plot.RData")
load(file="RData/L0_method3_venn_plot.RData")
load(file="RData/L1_venn_plot.RData")
load(file="RData/L1_method3_venn_plot.RData")
load(file="RData/MIO_venn_plot.RData")
load(file="RData/MIO_method3_venn_plot.RData")
load(file="RData/CLIME_venn_plot.RData")
CLIME_F_sigma <- F_omiga
load(file="RData/Glasso_venn_plot.RData")
Glasso_F_sigma <- F_omiga

#Dimension
p <- 300

#T(Omega_hat_S|Z0(Omega_hat_US),SL(Omega_hat_S)) for L0
L0 <- {}
#T(Omega_hat_S|Z0(T_hat),SL(Omega_hat_S)) for L1
L1 <- {}
#T(Omega_hat_S|Z0(Omega_hat_US),SL(Omega_hat_S)) for MIO
MIO <- {}
CLIME <- {}
Glasso <- {}
for(i in 2:p){
  for(j in 1:(i-1)){
    if(sdar_unbias_F_sigma_star[i,j]!=0){
      L0_k <- paste(rownames(sdar_unbias_F_sigma_star)[i],colnames(sdar_unbias_F_sigma_star)[j],sep="-")
      L0 <- c(L0,L0_k)
    }
  }
}

for(i in 2:p){
  for(j in 1:(i-1)){
    if(NS_T_star2[i,j]!=0){
      L1_k <- paste(rownames(NS_T_star2)[i],colnames(NS_T_star2)[j],sep="-")
      L1 <- c(L1,L1_k)
    }
  }
}

for(i in 2:p){
  for(j in 1:(i-1)){
    if(MIO_unbias_F_sigma_star[i,j]!=0){
      MIO_k <- paste(rownames(MIO_unbias_F_sigma_star)[i],colnames(MIO_unbias_F_sigma_star)[j],sep="-")
      MIO <- c(MIO,MIO_k)
    }
  }
}

for(i in 2:p){
  for(j in 1:(i-1)){
    if(CLIME_F_sigma[i,j]!=0){
      CLIME_k <- paste(rownames(CLIME_F_sigma)[i],colnames(CLIME_F_sigma)[j],sep="-")
      CLIME <- c(CLIME,CLIME_k)
    }
  }
}

for(i in 2:p){
  for(j in 1:(i-1)){
    if(Glasso_F_sigma[i,j]!=0){
      Glasso_k <- paste(rownames(Glasso_F_sigma)[i],colnames(Glasso_F_sigma)[j],sep="-")
      Glasso <- c(Glasso,Glasso_k)
    }
  }
}

x1 <- setdiff(L0,L1)
x2 <- setdiff(x1,CLIME)
x3 <- setdiff(x2,Glasso)
x4 <- setdiff(x3,MIO)
#46 direct connections only identified by T(Omega_hat_S|Z0(Omega_hat_US),SL(Omega_hat_S)) for L0
length(x4) #46
library(stringr)
probe_id <- str_split_fixed(x4,"-",n=2)
library(hgu133plus2.db)
ids <- toTable(hgu133plus2SYMBOL)
gene <- data.frame(matrix(NA,nrow=length(x4),ncol=5))
colnames(gene) <- c("gene1","gene2","partial","L1_Omega4_partial", "MIO_Omega4_partial")
gene$gene1 <- ids$symbol[match(probe_id[,1],ids$probe_id)]
gene$gene2 <- ids$symbol[match(probe_id[,2],ids$probe_id)]

#Estimated partial correlation of T(Omega_hat_S|Z0(Omega_hat_US),SL(Omega_hat_S)) for L0
sdar_partial_correlation <- sdar_unbias_F_sigma_star
for(i in 1:p){
  for(j in 1:p){
    if(i!=j){
      sdar_partial_correlation[i,j] <- -sdar_unbias_F_sigma_star[i,j]/(sqrt(sdar_unbias_F_sigma_star[i,i]*sdar_unbias_F_sigma_star[j,j]))
    }
    else{sdar_partial_correlation[i,j] <- 0}
  }
}
for(i in 1:length(x4)){
  gene$partial[i] <- sdar_partial_correlation[which(rownames(sdar_partial_correlation)==probe_id[i,1]),which(colnames(sdar_partial_correlation)==probe_id[i,2])]
}

#Estimated partial correlation for T(T_hat|Z0(T_hat),SL(T_hat)) for L1
NS_method3_partial_correlation <- NS_T_star1
for(i in 1:p){
  for(j in 1:p){
    if(i!=j){
      NS_method3_partial_correlation[i,j] <- -NS_T_star1[i,j]/(sqrt(NS_T_star1[i,i]*NS_T_star1[j,j]))
    }
    else{NS_method3_partial_correlation[i,j] <- 0}
  }
}
for(i in 1:length(x4)){
  gene$L1_Omega4_partial[i] <- NS_method3_partial_correlation[which(rownames(NS_method3_partial_correlation)==probe_id[i,1]),which(colnames(NS_method3_partial_correlation)==probe_id[i,2])]
}

#Estimated partial correlation for T(T_hat|Z0(T_hat),SL(T_hat)) for MIO
MIO_method3_partial_correlation <- MIO_T_star1
for(i in 1:p){
  for(j in 1:p){
    if(i!=j){
      MIO_method3_partial_correlation[i,j] <- -MIO_T_star1[i,j]/(sqrt(MIO_T_star1[i,i]*MIO_T_star1[j,j]))
    }
    else{MIO_method3_partial_correlation[i,j] <- 0}
  }
}
for(i in 1:length(x4)){
  gene$MIO_Omega4_partial[i] <- MIO_method3_partial_correlation[which(rownames(MIO_method3_partial_correlation)==probe_id[i,1]),which(colnames(MIO_method3_partial_correlation)==probe_id[i,2])]
}




gene
#15 in 46 which are not founded by T(T_hat|Z0(T_hat),SL(T_hat)) for L1 and MIO
length(which(gene$L1_Omega4_partial==0 & gene$MIO_Omega4_partial==0)) #15
gene[which(gene$L1_Omega4_partial==0 & gene$MIO_Omega4_partial==0),]
probe_id[which(gene$L1_Omega4_partial==0 & gene$MIO_Omega4_partial==0),]
gene_15 <- cbind(probe_id[which(gene$L1_Omega4_partial==0 & gene$MIO_Omega4_partial==0),],gene[which(gene$L1_Omega4_partial==0 & gene$MIO_Omega4_partial==0),])
gene_15_order <- gene_15[order(abs(gene_15$partial),decreasing=TRUE),]
gene_15_order

#Estimated absolute partial correlation matrix of T(Omega_hat_S|Z0(T_hat),SL(Omega_hat_S)) for L1
L1_partial_correlation <- NS_T_star2
for(i in 1:p){
  for(j in 1:p){
    if(i!=j){
      L1_partial_correlation[i,j] <- abs(-NS_T_star2[i,j]/(sqrt(NS_T_star2[i,i]*NS_T_star2[j,j])))
    }
    else{L1_partial_correlation[i,j] <- 0}
  }
}

#Estimated partial correlation of T(Omega_hat_S|Z0(Omega_hat_US),SL(Omega_hat_S)) for MIO
MIO_partial_correlation <- MIO_unbias_F_sigma_star
for(i in 1:p){
  for(j in 1:p){
    if(i!=j){
      MIO_partial_correlation[i,j] <- -MIO_unbias_F_sigma_star[i,j]/(sqrt(MIO_unbias_F_sigma_star[i,i]*MIO_unbias_F_sigma_star[j,j]))
    }
    else{MIO_partial_correlation[i,j] <- 0}
  }
}

#Estimated absolute partial correlation matrix of CLIME
CLIME_partial_correlation <- CLIME_F_sigma
for(i in 1:p){
  for(j in 1:p){
    if(i!=j){
      CLIME_partial_correlation[i,j] <- abs(-CLIME_F_sigma[i,j]/(sqrt(CLIME_F_sigma[i,i]*CLIME_F_sigma[j,j])))
    }
    else{CLIME_partial_correlation[i,j] <- 0}
  }
}

#Estimated absolute partial correlation matrix of GLasso
Glasso_partial_correlation <- Glasso_F_sigma
for(i in 1:p){
  for(j in 1:p){
    if(i!=j){
      Glasso_partial_correlation[i,j] <- abs(-Glasso_F_sigma[i,j]/(sqrt(Glasso_F_sigma[i,i]*Glasso_F_sigma[j,j])))
    }
    else{Glasso_partial_correlation[i,j] <- 0}
  }
}


x5 <- setdiff(L1,L0)
probe_L1 <- str_split_fixed(x5,"-",n=2)
probe_L1_partial <- {}
for(i in 1:nrow(probe_L1)){
  probe_L1_partial_k <- L1_partial_correlation[which(rownames(L1_partial_correlation)==probe_L1[i,1]),which(colnames(L1_partial_correlation)==probe_L1[i,2])]
  probe_L1_partial <- c(probe_L1_partial,probe_L1_partial_k)
}
gene_L1 <- data.frame(matrix(NA,nrow=length(x5),ncol=0))
gene_L1$probe <- x5
gene_L1$partial <- probe_L1_partial
gene_L1

x6 <- setdiff(CLIME,L0)
probe_CLIME <- str_split_fixed(x6,"-",n=2)
probe_CLIME_partial <- {}
for(i in 1:nrow(probe_CLIME)){
  probe_CLIME_partial_k <- CLIME_partial_correlation[which(rownames(CLIME_partial_correlation)==probe_CLIME[i,1]),which(colnames(CLIME_partial_correlation)==probe_CLIME[i,2])]
  probe_CLIME_partial <- c(probe_CLIME_partial,probe_CLIME_partial_k)
}
gene_CLIME <- data.frame(matrix(NA,nrow=length(x6),ncol=0))
gene_CLIME$probe <- x6
gene_CLIME$partial <- probe_CLIME_partial
gene_CLIME

x7 <- setdiff(Glasso,L0)
probe_Glasso <- str_split_fixed(x7,"-",n=2)
probe_Glasso_partial <- {}
for(i in 1:nrow(probe_Glasso)){
  probe_Glasso_partial_k <- Glasso_partial_correlation[which(rownames(Glasso_partial_correlation)==probe_Glasso[i,1]),which(colnames(Glasso_partial_correlation)==probe_Glasso[i,2])]
  probe_Glasso_partial <- c(probe_Glasso_partial,probe_Glasso_partial_k)
}
gene_Glasso <- data.frame(matrix(NA,nrow=length(x7),ncol=0))
gene_Glasso$probe <- x7
gene_Glasso$partial <- probe_Glasso_partial
gene_Glasso

x8 <- setdiff(MIO,L0)
probe_MIO <- str_split_fixed(x8,"-",n=2)
probe_MIO_partial <- {}
for(i in 1:nrow(probe_MIO)){
  probe_MIO_partial_k <- MIO_partial_correlation[which(rownames(MIO_partial_correlation)==probe_MIO[i,1]),which(colnames(MIO_partial_correlation)==probe_MIO[i,2])]
  probe_MIO_partial <- c(probe_MIO_partial,probe_MIO_partial_k)
}
gene_MIO <- data.frame(matrix(NA,nrow=length(x8),ncol=0))
gene_MIO$probe <- x8
gene_MIO$partial <- probe_MIO_partial
gene_MIO



probe_total <- unique(c(x5,x6,x7,x8))
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
  index <- which(gene_MIO$probe==gene_total$probe[i])
  if(length(index)==0){
    gene_total$MIO[i] <- 0
  }
  else{gene_total$MIO[i] <- gene_MIO$partial[index]}
}


for(i in 1:nrow(gene_total)){
  value_vector <- c(gene_total$L1[i],gene_total$CLIME[i],gene_total$Glasso[i],gene_total$MIO[i])
  gene_total$secondsmall[i] <- value_vector[order(value_vector)][3]
}

gene_total
#18 have absolute partial correlation estimate larger than 0.1 by at least two other methods
length(which(gene_total$secondsmall>0.1)) #18

#Estimated absolute partial correlation matrix of T(T_hat|Z0(T_hat),SL(T_hat)) for L0
L0_omega4_partial_correlation <- sdar_T_star1
for(i in 1:p){
  for(j in 1:p){
    if(i!=j){
      L0_omega4_partial_correlation[i,j] <- abs(-sdar_T_star1[i,j]/(sqrt(sdar_T_star1[i,i]*sdar_T_star1[j,j])))
    }
    else{L0_omega4_partial_correlation[i,j] <- 0}
  }
}

gene_18 <- gene_total[which(gene_total$secondsmall>0.1),]
gene_18_probesplit <- str_split_fixed(gene_18$probe,"-",n=2) 
for(i in 1:nrow(gene_18)){
  gene_18$L0_omega4_partial[i] <- L0_omega4_partial_correlation[which(rownames(L0_omega4_partial_correlation)==gene_18_probesplit[i,1]),which(colnames(L0_omega4_partial_correlation)==gene_18_probesplit[i,2])]
}
#Only 3 of 18 is detected by T(T_hat|Z0(T_hat),SL(T_hat)) for L0
length(which(gene_18$L0_omega4_partial!=0)) #3


#Venn diagram
x <- list(
  A = L0,
  B = L1, 
  C = CLIME,
  D = Glasso,
  E = MIO
)

display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, cat.dist=c(0.20,0.20,0.24,0.21,0.20), cat.pos = c(50, 320, 250, 150, 35), cat.cex=1.2, ...)
  grid.draw(venn_object)
}

#Venn diagram for estimated absolute partial correlations
display_venn(
  x,
  category.names = c(expression(L[0]) , expression(L[1]) , "CLIME", "GLasso", "MIO"),
  fill = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#FF6633")
)


load(file="RData/L0_venn_plot_threshold1.RData")
load(file="RData/L1_venn_plot_threshold1.RData")
load(file="RData/CLIME_venn_plot_threshold1.RData")
load(file="RData/Glasso_venn_plot_threshold1.RData")
load(file="RData/MIO_venn_plot_threshold1.RData")

#Dimension
p <- 300

L0_1 <- {}
L1_1 <- {}
CLIME_1 <- {}
Glasso_1 <- {}
MIO_1 <- {}
for(i in 2:p){
  for(j in 1:(i-1)){
    if(sdar_partial_correlation_threshold1[i,j]!=0){
      L0_1k <- paste(rownames(sdar_partial_correlation_threshold1)[i],colnames(sdar_partial_correlation_threshold1)[j],sep="-")
      L0_1 <- c(L0_1,L0_1k)
    }
  }
}

for(i in 2:p){
  for(j in 1:(i-1)){
    if(NS_partial_correlation_threshold1[i,j]!=0){
      L1_1k <- paste(rownames(NS_partial_correlation_threshold1)[i],colnames(NS_partial_correlation_threshold1)[j],sep="-")
      L1_1 <- c(L1_1,L1_1k)
    }
  }
}

for(i in 2:p){
  for(j in 1:(i-1)){
    if(MIO_partial_correlation_threshold1[i,j]!=0){
      MIO_1k <- paste(rownames(MIO_partial_correlation_threshold1)[i],colnames(MIO_partial_correlation_threshold1)[j],sep="-")
      MIO_1 <- c(MIO_1,MIO_1k)
    }
  }
}

for(i in 2:p){
  for(j in 1:(i-1)){
    if(CLIME_partial_correlation_threshold1[i,j]!=0){
      CLIME_1k <- paste(rownames(CLIME_partial_correlation_threshold1)[i],colnames(CLIME_partial_correlation_threshold1)[j],sep="-")
      CLIME_1 <- c(CLIME_1,CLIME_1k)
    }
  }
}

for(i in 2:p){
  for(j in 1:(i-1)){
    if(Glasso_partial_correlation_threshold1[i,j]!=0){
      Glasso_1k <- paste(rownames(Glasso_partial_correlation_threshold1)[i],colnames(Glasso_partial_correlation_threshold1)[j],sep="-")
      Glasso_1 <- c(Glasso_1,Glasso_1k)
    }
  }
}

x <- list(
  A = L0_1,
  B = L1_1, 
  C = CLIME_1,
  D = Glasso_1,
  E = MIO_1
)

display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, cat.dist=c(0.20,0.20,0.24,0.21,0.20), cat.pos = c(50, 320, 250, 150, 35), cat.cex=1.2, ...)
  grid.draw(venn_object)
}

#Venn diagram for estimated absolute partial correlations thresholded by 0.1
display_venn(
  x,
  category.names = c(expression(L[0]) , expression(L[1]) , "CLIME", "GLasso", "MIO"),
  fill = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#FF6633")
)




