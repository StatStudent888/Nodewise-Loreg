library("limma")
#Load data
MDA<-read.csv("MDA133.csv",sep = ",")
MDA_L<-read.csv("MDA133_L.csv",sep = ",")


#View missing values
sum(is.na(MDA))
summary(MDA[,-1])


gene_names<-MDA[,1]  
sample_names<-colnames(MDA)[-1] 
gene_MDA<-MDA[,-1]
rownames(gene_MDA)=gene_names
colnames(gene_MDA)=sample_names

rownames(MDA_L)<-sample_names
MDA_L_PCR<-MDA_L[which(MDA_L[,3]==1),]
PCR_name<-rownames(MDA_L_PCR)
MDA_L_RD<-MDA_L[which(MDA_L[,3]==0),]
RD_name<-rownames(MDA_L_RD)

#Log transformation
MDA_log2<-log2(gene_MDA+1)
summary(MDA_log2)


################################
#Boxplot
################################
#Before log transformation
par(cex = 0.3)  
n.sample=ncol(gene_MDA)
cols <- rainbow(n.sample*1.2) 
dev.off()
boxplot(gene_MDA, col = cols,main=" gene_MDA expression value",las=2)

#After log transformation
par(cex = 0.3)  
n.sample=ncol(MDA_log2)
cols <- rainbow(n.sample*1.2) 
dev.off()
boxplot(MDA_log2, col = cols,main="MDA_log2 expression value",las=2)


################################
#Screening for differentially expressed genes
################################


###################
#   1    Construct design matrix
####################

group_list=as.character(MDA_L[,2])
suppressMessages(library(limma)) 
design <- model.matrix(~0+factor(group_list))   
colnames(design)=levels(factor(group_list))                                             
rownames(design)=sample_names             
design 

#####################################
#   2    Construct contrast matrix
#####################################
contrast.matrix<-makeContrasts(paste0(unique(group_list),collapse = "-"),levels = design)
contrast.matrix

#####################################
#   3     Linear model analysis
#####################################
fit <- lmFit(MDA_log2,design) 
fit1 <- contrasts.fit(fit, contrast.matrix) 
fit2 <- eBayes(fit1)
Output1 = topTable(fit2, coef=1,sort.by="P", n=Inf)  
sum(is.na(Output1))
nrDEG1 = na.omit(Output1)

#Select 300 differentially expressed genes
Test_log2_genename<-rownames(nrDEG1[1:300,])
Test_log2_genedata<-MDA_log2[Test_log2_genename,]
Test_log2_genedata<-as.data.frame(Test_log2_genedata)
PCR_Test_log2_genedata<-Test_log2_genedata[,PCR_name]
RD_Test_log2_genedata<-Test_log2_genedata[,RD_name]


write.csv(PCR_Test_log2_genedata,file = "Log2/PCR_Test_log2_genedata_300.csv")
write.csv(RD_Test_log2_genedata,file = "Log2/RD_Test_log2_genedata_300.csv")
