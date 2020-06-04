
rm(list=ls())
setwd("D:/1A_MM诊断模型构建_20200510")


# Read the counts of 37 samples (kraken，non-standardization)
Species_train <- read.csv("Species_train.csv", header=T, row.names=1)
# Species_test <- read.csv("Species_test.csv", header=T, row.names=1)


species_DNA <- Species_train #######
species_DNA_HC <- species_DNA[,1:18]
species_DNA_MM <- species_DNA[,19:37]

### 1. core microbiome: >0.0001 = 0.01% = 10e-4
noise.removal <- function(dataframe, percent=0.01, top=NULL){
  dataframe->Matrix
  bigones <- rowSums(Matrix)/(sum(rowSums(Matrix))) >= percent
  Matrix_1 <- Matrix[bigones,]
  print(percent)
  return(Matrix_1)}

species_DNA_HC_1 <- noise.removal(species_DNA_HC, percent=0.0001) # 123
species_DNA_MM_1 <- noise.removal(species_DNA_MM, percent=0.0001) # 153

### 2. core microbiome: >=90% samples
nonzero.sample <- function(dataframe, detection=0.9, top=NULL){
  dataframe -> Matrix
  index = c()
  for (i in 1:nrow(Matrix)){
    if (length(which(Matrix[i,] != 0))/length(Matrix[i,]) >= detection)  #
      index = c(index,i)}
  Matrix[index,] -> Matrix_1
}

species_DNA_HC_2 <- nonzero.sample(species_DNA_HC_1, detection=0.9) # 123
species_DNA_MM_2 <- nonzero.sample(species_DNA_MM_1, detection=0.9) # 153

write.csv(species_DNA_HC_2, "Figures/core_microbiome_species_HC.csv")
write.csv(species_DNA_MM_2, "Figures/core_microbiome_species_MM.csv")



## core microbiome, Venn diagram
# install.packages("VennDiagram")
library(VennDiagram)
venn.diagram(list(HC=row.names(species_DNA_HC_2),
                  MM=row.names(species_DNA_MM_2)),
             cex=2, fontfamily=1,
             cat.cex=2.5, cat.fontface=1, 
             fill=c("DodgerBlue","red"), alpha=c(0.5,0.5),scaled=TRUE,
             filename="Figures/core_microbiome_species.tiff",
             main="Core bacterial species", main.cex=2)



# DESeq2, identification of differential species
library(DESeq2) #
a <- union(row.names(species_DNA_HC_2), row.names(species_DNA_MM_2)) # 184
Abundance_S1 <- Species_train[a,]
data <- Abundance_S1 # count reads
group <- as.factor(c(rep("HC",18),rep("MM",19)))
conditions=data.frame(colnames(data),group)
dds <- DESeqDataSetFromMatrix(data, colData=conditions, design= ~group)
dds <- dds[ rowSums(counts(dds)) > 10,]
dds2 <- DESeq(dds) #
DATA_normalization <- as.data.frame(counts(dds2, normalized=T))
# write.csv(DATA_normalization, "Differential species/DATA_species_DESeq2.csv")

res<-results(dds2, contrast=c("group","MM","HC")) 
deg2=as.data.frame(res)
deg2 <- deg2[order(deg2$log2FoldChange,decreasing=T), ]
resdata <- merge(deg2, DATA_normalization, by="row.names", sort=F)
result_DESeq2=resdata[which(resdata$padj< 0.01 & abs(resdata$log2FoldChange)> 1.5),] # 177
result_DESeq2=result_DESeq2[order(result_DESeq2$log2FoldChange,decreasing=T),]
write.csv(result_DESeq2,"Figures/Diff_species_DESeq2.csv",row.names=F) # 



#########################################################################
# 38 species =10 MM-depleted + 28 MM-enriched
# visualization
result_DESeq2 <- read.csv("Figures/Diff_species_DESeq2_selected.csv",header=T,row.names=1)
Diff_species <- result_DESeq2[,1:6]
Diff_species_abundance <- result_DESeq2[,-c(1:6)]
DATA_0 <- log10(Diff_species_abundance + 1)
DATA <- t(scale(t(DATA_0), center=T, scale=T))


# OTU <- read.csv("DNA_Bacteria.csv",header=T,row.names=1)
# A <- paste("s__", Diff_species[,1], sep="")
# Num <- c()
# for (i in 1:length(A)) {
#   num <- grep(A[i], row.names(OTU))
#   Num <- c(Num, num[1])
# }
# OTU_s <- OTU[Num,]
# write.csv(OTU_s, "OTU_s.csv")

OTU_s <- read.csv("OTU_s.csv", header=F, row.names=7)
OTU_s <- OTU_s[row.names(DATA),]

annotation_col = data.frame(Group=as.factor(c(rep("HC",18),rep("MM",19))))
rownames(annotation_col) = colnames(DATA)
annotation_row = data.frame(Family=as.factor(OTU_s$V5), Phylum=as.factor(OTU_s$V2))
rownames(annotation_row) = rownames(DATA)

ann_colors = list(Group=c(HC="blue2", MM="red2"),
                  Phylum=c("p__Actinobacteria"="Purple", "p__Bacteroidetes"="SkyBlue",
                           "p__Fusobacteria"="Gray", "p__Firmicutes"="Cyan",
                           "p__Proteobacteria"="Orange"))

c1=colorRampPalette(c("blue","gray"))(-min(DATA)*1000)
c2=colorRampPalette(c("gray","red"))(max(DATA)*1000)

library(pheatmap)
pheatmap(DATA,cellwidth=12,cellheight=12,border_color=NA,cluster_rows=T,cluster_col=F,
         scale="none",show_rownames=T,show_colnames=F,color=c(c1,c2),angle_col=315, 
         fontsize=10,fontsize_row=10,fontsize_col=10,fontsize_number=10, 
         clustering_distance_rows="correlation",clustering_distance_cols="correlation",
         annotation_col=annotation_col,annotation_row=annotation_row,
         annotation_colors=ann_colors,
         filename="Figures/Diff_species_heatmap_selected_2.pdf",width=15,height=12)

# pheatmap(DATA,cellwidth=12,cellheight=12,border_color=NA,cluster_rows=T,cluster_col=F,
#          scale="none",show_rownames=T,show_colnames=F,color=c(c1,c2),angle_col=315, 
#          fontsize=10,fontsize_row=10,fontsize_col=10,fontsize_number=10, 
#          clustering_distance_rows="manhattan",clustering_distance_cols="correlation",
#          annotation_col=annotation_col,annotation_row=annotation_row,
#          annotation_colors=ann_colors,
#          filename="Figures/Diff_species_heatmap_selected.pdf",width=15,height=12)




