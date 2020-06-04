
rm(list=ls())
setwd("D:/1A_MM诊断模型构建_20200510")

library(caret)
library(nnet) # 神经网络 
library(randomForest) # 随机森林 rf
library(arm) # 贝叶斯广义线性模型 
library(pls) # 偏最小二乘法 pls
library(kknn) # knn
library(kernlab) # SVM
library(pamr) # pam
library(klaR) # rda
library(ROCR)
library(pROC)
library(fastAdaboost) # adaboost


# 11 species were combineded as features：
Features<-c("Klebsiella_pneumoniae","Klebsiella_variicola","Klebsiella_aerogenes",
            "Enterobacter_cloacae","Citrobacter_freundii",
            "Streptococcus_oralis","Streptococcus_gordonii","Streptococcus_salivarius",
            "Streptococcus_mitis","Clostridium_butyricum","Anaerostipes_hadrus") 

### training dataset
Species_train <- read.csv("Species_train.csv", header=T, row.names=1)
Species_train <- Species_train[Features,]
#Data_train <- data.frame(t(Species_train), group=as.factor(c(rep("HC",18), rep("MM",19))))
trainx <- t(Species_train)
trainy <- as.factor(c(rep("HC",18), rep("MM",19)))
write.csv(trainx, "trainx.csv")
trainx <- read.csv("trainx.csv", header=T, row.names=1)


### validation dataset
Species_test <- read.csv("Species_test.csv", header=T, row.names=1)
Species_test <- Species_test[Features,]
#Data_test <- data.frame(t(Species_test), group=as.factor(c(rep("HC",3), rep("MM",11))))
testx <- t(Species_test)
testy <- as.factor(c(rep("HC",3), rep("MM",11)))
write.csv(testx, "testx.csv")
testx <- read.csv("testx.csv", header=T, row.names=1)

# training models
set.seed(123)
fitControl <- trainControl(method="repeatedcv",number=5,repeats=3,classProbs=T,
                           summaryFunction=twoClassSummary,savePredictions=T)

modelfit1 <- train(trainx, trainy, method='rf', trControl=fitControl, metric="ROC", verbose=F)
modelfit2 <- train(trainx, trainy, method="adaboost", trControl=fitControl, metric="ROC", verbose=F)
modelfit3 <- train(trainx, trainy, method="nnet", trControl=fitControl, metric="ROC", verbose=F)
modelfit4 <- train(trainx, trainy, method="simpls", trControl=fitControl, metric="ROC", verbose=F)
modelfit5 <- train(trainx, trainy, method="kknn", trControl=fitControl, metric="ROC")
modelfit6 <- train(trainx, trainy, method="bayesglm", trControl=fitControl, metric="ROC")
modelfit7 <- train(trainx, trainy, method="lda", trControl=fitControl, metric="ROC", verbose=F)
modelfit8 <- train(trainx, trainy, method="pam", trControl=fitControl, metric="ROC")
modelfit9 <- train(trainx, trainy, method="svmLinear",trControl=fitControl,tuneLength=8,metric="ROC",verbose=F)
modelfit10 <- train(trainx, trainy, method="svmRadial",trControl=fitControl,tuneLength=8,metric="ROC",verbose=F)


###  10 fitting models, their characteristics：ROC, Sens, Spec
resamps <- resamples(list(modelfit1, modelfit2, modelfit3, modelfit4, modelfit5,
                         modelfit6, modelfit7, modelfit8, modelfit9, modelfit10))

#summary(resamps)
#splom(resamps)


pdf("Figures/Characteristics_of_models.pdf",width=6, height=5)
theme1 <- trellis.par.get()
theme1$plot.symbol$col = c()
theme1$plot.symbol$pch = 1
theme1$plot.line$col = c()
theme1$plot.line$lwd <- 2
trellis.par.set(theme1)
bwplot(resamps, layout=c(3, 1), box.ratio=2, colors="red")
dev.off()


pdf("Figures/ROC_of_models.pdf",width=6, height=6)
trellis.par.set(caretTheme())
dotplot(resamps, metric="ROC", box.ratio=2)
dev.off()




# Calculates cross-tabulation of observed and predicted classes
models <- list(modelfit1, modelfit2, modelfit3, modelfit4, modelfit5,
               modelfit6, modelfit7, modelfit8, modelfit9, modelfit10)

predValues = extractPrediction(models, testX=testx, testY=testy)
testValues = subset(predValues, dataType=="Test")
Pred1 = subset(testValues, model=="rf")
Pred2 = subset(testValues, model=="svmRadial")
Pred3 = subset(testValues, model=="nnet")
Pred4 = subset(testValues, model=="simpls")
Pred5 = subset(testValues, model=="kknn")
Pred6 = subset(testValues, model=="bayesglm")
Pred7 = subset(testValues, model=="lda")
Pred8 = subset(testValues, model=="rda")
Pred9 = subset(testValues, model=="svmLinear")
Pred10 = subset(testValues, model=="adaboost")
confusionMatrix(Pred1$pred, Pred1$obs)
confusionMatrix(Pred2$pred, Pred2$obs)
confusionMatrix(Pred3$pred, Pred3$obs)
confusionMatrix(Pred4$pred, Pred4$obs)
confusionMatrix(Pred5$pred, Pred5$obs)
confusionMatrix(Pred6$pred, Pred6$obs)
confusionMatrix(Pred7$pred, Pred7$obs)
confusionMatrix(Pred8$pred, Pred8$obs)
confusionMatrix(Pred9$pred, Pred9$obs)
confusionMatrix(Pred10$pred, Pred10$obs)


### plot
probValues = extractProb(models, testX=testx, testY=testy)
testProbs = subset(probValues, dataType=="Test")
prob1 = subset(testProbs, model=="rf")
prob2 = subset(testProbs, model=="adaboost")
prob3 = subset(testProbs, model=="nnet")
prob4 = subset(testProbs, model=="simpls")
prob5 = subset(testProbs, model=="kknn")
prob6 = subset(testProbs, model=="bayesglm")
prob7 = subset(testProbs, model=="lda")
prob8 = subset(testProbs, model=="pam")
prob9 = subset(testProbs, model=="svmLinear")
prob10 = subset(testProbs, model=="svmRadial")

prob1$lable=ifelse(prob1$obs=='MM', yes=1, 0)
prob2$lable=ifelse(prob2$obs=='MM', yes=1, 0)
prob3$lable=ifelse(prob3$obs=='MM', yes=1, 0)
prob4$lable=ifelse(prob4$obs=='MM', yes=1, 0)
prob5$lable=ifelse(prob5$obs=='MM', yes=1, 0)
prob6$lable=ifelse(prob6$obs=='MM', yes=1, 0)
prob7$lable=ifelse(prob7$obs=='MM', yes=1, 0)
prob8$lable=ifelse(prob8$obs=='MM', yes=1, 0)
prob9$lable=ifelse(prob9$obs=='MM', yes=1, 0)
prob10$lable=ifelse(prob10$obs=='MM', yes=1, 0)

pred1 = prediction(prob1$MM, prob1$lable)
pred2 = prediction(prob2$MM, prob2$lable)
pred3 = prediction(prob3$MM, prob3$lable)
pred4 = prediction(prob4$MM, prob4$lable)
pred5 = prediction(prob5$MM, prob5$lable)
pred6 = prediction(prob6$MM, prob6$lable)
pred7 = prediction(prob7$MM, prob7$lable)
pred8 = prediction(prob8$MM, prob8$lable)
pred9 = prediction(prob9$MM, prob9$lable)
pred10 = prediction(prob10$MM, prob10$lable)

perf1 = performance(pred1, measure="tpr", x.measure="fpr")
perf2 = performance(pred2, measure="tpr", x.measure="fpr")
perf3 = performance(pred3, measure="tpr", x.measure="fpr")
perf4 = performance(pred4, measure="tpr", x.measure="fpr")
perf5 = performance(pred5, measure="tpr", x.measure="fpr")
perf6 = performance(pred6, measure="tpr", x.measure="fpr")
perf7 = performance(pred7, measure="tpr", x.measure="fpr")
perf8 = performance(pred8, measure="tpr", x.measure="fpr")
perf9 = performance(pred9, measure="tpr", x.measure="fpr")
perf10 = performance(pred10, measure="tpr", x.measure="fpr")




color <- rainbow(10)
methods <- c("rf","adaboost","nnet","simpls","kknn","bayesglm","lda","pam","svmLines","svmRadial")



i=10
pdf(paste0("Test_pridict_models/model_",i,".pdf",sep=""),width=6, height=6)

plot(perf10, col=color[1], main=paste("Model0",i,sep=""),
     xlab='False positive rate', ylab='True positive rate',
     lwd=2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)

abline(0, 1, lty=2)
legend("bottomright", col=color[1], lty=1, lwd=2,
       legend=c(paste(methods[i],' AUC = ',round(as.numeric(performance(pred10,'auc')@y.values),3))))
                
dev.off()








# #tiff("Test_pridict_models_1.tiff", width=15, height=15, units="cm", res=300)
# pdf("Test_pridict_models_1.pdf",width=6, height=6)
# plot(perf1, col='gray', xlab='False positive rate', ylab='True positive rate',
#      lwd=2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2, add=F)
# plot(perf2, col='orange', lwd=2, add=T)
# plot(perf3, col='green', lwd=2, add=T)
# plot(perf4, col='blue', lwd=2, add=T)
# plot(perf5, col='red', lwd=2, add=T)
# 
# abline(0,1,lty=2) 
# legend("bottomright",
#        legend=c(paste('rf AUC = ',round(as.numeric(performance(pred1,'auc')@y.values),3)),
#                 paste('adaboost AUC = ',round(as.numeric(performance(pred2,'auc')@y.values),3)),
#                 paste('nnet AUC = ',round(as.numeric(performance(pred3,'auc')@y.values),3)),
#                 paste('simpls AUC = ',round(as.numeric(performance(pred4,'auc')@y.values),3)),
#                 paste('kknn AUC = ',round(as.numeric(performance(pred5,'auc')@y.values),3))),
#        col=c('gray','orange','green','blue','red'),lty=1,lwd=2)
# dev.off()
# 
# 
# #tiff("Test_pridict_models_2.tiff", width=15, height=15,units="cm", res=300)
# pdf("Test_pridict_models_2.pdf",width=6, height=6)
# plot(perf6, col='gray', xlab='False positive rate', ylab='True positive rate',
#      lwd=2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2, add=F)
# plot(perf7, col='orange', lwd=2, add=T)
# plot(perf8, col='green', lwd=2, add=T)
# plot(perf9, col='blue', lwd=2, add=T)
# plot(perf10, col='red', lwd=2, add=T)
# 
# abline(0,1,lty=2) 
# legend("bottomright",
#        legend=c(paste('bayesglm AUC = ',round(as.numeric(performance(pred6,'auc')@y.values),3)),
#                 paste('lda AUC = ',round(as.numeric(performance(pred7,'auc')@y.values),3)),
#                 paste('pam AUC = ',round(as.numeric(performance(pred8,'auc')@y.values),3)),
#                 paste('svmLines AUC = ',round(as.numeric(performance(pred9,'auc')@y.values),3)),
#                 paste('svmRadial AUC = ',round(as.numeric(performance(pred10,'auc')@y.values),3))),
#        col=c('gray','orange','green','blue','red'),lty=1,lwd=2)
# dev.off()



