#random Forests fun
#Written by Jason Li

library(randomForest)

training.highCV<-trainExpressionData.mod.reorder[highCV,] #subset of high CV genes

trainExprsub.pca.reorder<-trainExprsub.pca[testedCellLines.mod,]

bestmtry <- tuneRF(t(as.data.frame(trainExpressionData.mod.reorder)),as.factor(trainKey$Paclitaxel), ntreeTry=100, 
                   stepFactor=1.5,improve=0.01, trace=TRUE, plot=TRUE, dobest=FALSE)
#tuning RF

rf1<-randomForest(formula, trainingData, ntree=10501, importance=TRUE, keep.forest=TRUE)
#this one does a regression because of the formula

j=2
nclass<-sum(trainKey[,j])
priorclass<-sum(trainKey[,j])/25
library(base)
for (j in 2:13){
  rf2<-randomForest(t(as.data.frame(training.highCV)), 
                    y=as.factor(trainKey[,j]), 
                    mtry=40, #nodesize=8,
                    ntree=5501, keep.forest=TRUE, 
                    importance=TRUE)
  print(rf2$confusion)
}
#classwt=c(1-priorclass, priorclass))
#RF2 does classification
#using high CV genes
#change parameters mtry, nodesize, ntree

for (j in 2:13){
  rf3<-randomForest(as.data.frame(t(train.pca$x)), 
                    y=as.factor(trainKey[,j]), 
                    ntree=1501, keep.forest=TRUE, 
                    importance=TRUE) 
  print(rf3$confusion)
}
#RF3 is classification based on PCA

#Code to implement ROC evaluation
#adult.rf.pr = predict(adult.rf,type="prob",newdata=data$val)[,2]
#adult.rf.pred = prediction(adult.rf.pr, data$val$income)
#adult.rf.perf = performance(adult.rf.pred,"tpr","fpr")
#plot(adult.rf.perf,main="ROC Curve for Random Forest",col=2,lwd=2)
#abline(a=0,b=1,lwd=2,lty=2,col="gray")


#importance(rf2) messy output!
varImpPlot(rf2) #graphs Mean Decrease in Accuracy and Mean GINI
varImpPlot(rf3)

MDSplot(rf2, trainKey[,j])