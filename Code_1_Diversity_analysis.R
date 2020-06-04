
rm(list=ls())
setwd("D:/1A_MM诊断模型构建_20200503")

# Read the counts of 37 samples (kraken，non-standardization)
Species_train <- read.csv("Species_train.csv", header=T, row.names=1)
# Species_test <- read.csv("Species_test.csv", header=T, row.names=1)


##################################################################################
# two groups (HC-18, MM-19)
## beta diversity
library(ade4)
library(vegan)

data = t(Species_train) 
gp=as.factor(c(rep("HC",18), rep("MM",19)))

### Compute species dissimilarity between each pair of plots
dis <- vegdist(wisconsin(sqrt(data)), method="bray")
adonis2(dis~gp,gp,permutations=999,method="bray") # 0.001 **  (PERMANOVA)
res_anosim_bray <- anosim(dis, gp, permutations=999, distance="bray")
summary(res_anosim_bray) # 0.001, Analysis of similarities (ANOSIM)

pdf("Figures/Bray-Curtis Dissmilarity1.pdf", width=5, height=6)
plot(res_anosim_bray, col=c("gray","blue","red"), cex.axis=1.3, cex.lab=1.3,
     xlab="", ylab="Bray-Curtis Dissmilarity")
dev.off()

pdf("Figures/Distribution of Bray-Curtis dissimilarity111.pdf", width=6, height=6)
res_betadisper_bray<-betadisper(dis, group=gp)
plot(res_betadisper_bray, ellipse=T, hull=F, xlim=c(-0.3, 0.2), ylim=c(-0.2, 0.2),
     cex=1.2,label.cex=1,col=c("HC"='blue', "MM"='red'), pch=c("HC"=16, "MM"=17))
dev.off()


##############################################################################
# two groups (HC-18, MM-19)
# a-diversity
library(vegan)
data = t(Species_train) #
H <- diversity(data)
shannon <- diversity(data, "shannon")
simp <- diversity(data, "simpson")
invsimp <- diversity(data, "inv")
## Unbiased Simpson (Hurlbert 1971, eq. 5) with rarefy:
unbias.simp <- rarefy(data, 2) - 1
## Fisher alpha
alpha <- fisher.alpha(data)
## data richness (S) and Pielou's evenness (J):
S <- specnumber(data) ## rowSums(BCI > 0) does the same...
J <- H/log(S)

AAA <- cbind.data.frame(shannon,simp,invsimp,unbias.simp,alpha,J)
BBB <- t(estimateR(data))
Diversity <- cbind.data.frame(AAA, BBB[,c(1,2,4)])

colnames(Diversity)<-c("Shannon","Simpson","InvSimpson","unbias.Simpson","Fisher's_alpha",
                       "Pielou's_evenness","S.observed","S.chao1","S.ACE")
Diversity$Class <- c(rep("HC",18),rep("MM",19))
write.csv(Diversity,"Figures/SampleDiversity_species.csv")



result_HC=matrix(0,nrow=9,ncol=6)
result_MM=matrix(0,nrow=9,ncol=6)
rownames(result_HC) <- colnames(Diversity)[1:9]
colnames(result_HC) <- c("Min","1st Qu.","Median","Mean","3rd Qu.","Max") # quantile分位数
row.names(result_MM) <- colnames(Diversity)[1:9]
colnames(result_MM) <- c("Min","1st Qu.","Median","Mean","3rd Qu.","Max") # quantile分位数
for (i in 1:9){
  result_HC[i,1:6] <- summary(Diversity[1:18,i])
  result_MM[i,1:6] <- summary(Diversity[19:37,i])
}
result <- cbind.data.frame(result_HC, result_MM)


### alpha-diversity, species richness
Diversity_1 <- Diversity
Result <- NULL
for (i in 1:9) {
  wil <- data.frame(value=as.numeric(Diversity_1[,i]), class=factor(Diversity_1$Class))
  wil.test <- wilcox.test(value ~ class, data=wil, alternative="two.sided", exact=F, correct=F)
  t.test <- t.test(value ~ class, data=wil, alternative="two.sided", exact=F, correct=F)
  Result$HC_mean[i] <- mean(Diversity[1:18,i])
  Result$HC_sd[i] <- sd(Diversity[1:18,i])
  Result$MM_mean[i] <- mean(Diversity[19:37,i])
  Result$MM_sd[i] <- sd(Diversity[19:37,i])
  Result$wil.Pvalue[i] <- wil.test$p.value
  Result$t.Pvalue[i] <- t.test$p.value
}
Result <- data.frame(Result)
RESULT <- cbind.data.frame(result, Result)
write.csv(RESULT,"Figures/Alpha_diversity_species.csv")





