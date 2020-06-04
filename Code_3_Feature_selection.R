
rm(list=ls())
setwd("D:/1A_MM诊断模型构建_20200510")


# Read the counts of 37 samples (kraken，non-standardization)
Species_train <- read.csv("Species_train.csv", header=T, row.names=1)
Species_test <- read.csv("Species_test.csv", header=T, row.names=1)


# 11 species were considered as candidate feature：
Features<-c("Klebsiella_pneumoniae","Klebsiella_variicola","Klebsiella_aerogenes",
            "Enterobacter_cloacae","Citrobacter_freundii",
            "Streptococcus_oralis","Streptococcus_gordonii","Streptococcus_salivarius",
            "Streptococcus_mitis","Clostridium_butyricum","Anaerostipes_hadrus")  

Features <- as.character(Features) #
data <- t(Species_train[Features, ])
group <- as.factor(c(rep("HC",18), rep("MM",19)))
Data <- data.frame(data, group)


#################################################################################################
# Rank Features By Importance
library("caret")
set.seed(12) # 800  caretFuncs
# prepare training scheme
control <- trainControl(method="repeatedcv", number=5, repeats=3)
# train the model
model <- train(group~., data=Data, method="lvq", preProcess="scale", trControl=control)
# estimate variable importance
importance <- varImp(model, scale=FALSE)
# summarize importance
print(importance)

pdf("Figures/Features_importance.pdf",width=5, height=5)
plot(importance, xlim=c(0.5,1))
dev.off()




