#PCA
train.pca<-prcomp(t(trainExpressionData.mod.reorder))
#train.pca$x
#train.pca$rotation
summary(train.pca)

#PCA plot colored by response drug
# 2 is first drug CGC.11047
#black is no response, red is yes response
for (j in 2:13){
plot(train.pca$x[,1], train.pca$x[,2], col=(as.integer(trainKey[,j]+1)), pch=20, xlab="PCA1", ylab="PCA2", main=paste("Response to", names(trainKey)[j]))
text(train.pca$x[,1], train.pca$x[,2], rownames(train.pca$x), col=as.integer(trainKey[,j]+1))
}

# SVM Parameters 
svmCost     <- 20
svmGamma    <- 0.078
svmKernel   <- "linear"
svmDegree   <- 2
svmType     <- "C-classification"
svmCoef0    <- 1
svmCross    <- 1

# Select on High Variance Genes
SelectCV    <- TRUE
ThreshCV    <- 28
AbsValCV    <- TRUE
AboveThresh <- TRUE

# Select on Maximum Expressed Genes
SelectMax   <- FALSE
MaxThresh   <- 2
AbsValMax   <- TRUE

#Select on PCA
SelectPCA <-FALSE
# If SelectCV == FALSE, Use nGenes Selection Criteria
nStart      <- 1
nGenes      <- 1000

# General Parameters
drugs               <- gsub("-", ".", trainKeyTransposed$Drug) # Vector of Drug Names
nDrugs              <- length(drugs)
totalGeneCount      <- dim(as.data.frame(trainExpressionData))[1] # This is the Maximum Number of Predictors Possible

if (SelectCV == TRUE) {
  highCV          <- getByCV(trainingData, ThreshCV, 13, AbsValCV, AboveThresh)
  nGenes          <- length(highCV)
  predictorGenes  <- paste(paste("X", highCV, sep=""), collapse= " + ") #Creates Predictor String for Formula 
} else if (SelectMax == TRUE) {
  highExpression  <- getByCV(trainingData, MaxThresh, 13, AbsValMax)
  nGenes          <- length(highExpression)
  predictorGenes  <- paste(paste("X", highExpression, sep=""), collapse= " + ") #Creates Predictor String for Formula 
} else if (SelectPCA == TRUE) {
  nGenes <- 2
  predictorGenes <- ("PCA1+PCA2")
} else {
  genePredictorRange  <- nStart:(nStart + nGenes - 1) #Selects Genes Labeled X1 - X10 in Training Data Set - 100 genes, Trying not to Get Too Crazy
  predictorGenes  <- paste(paste("X", genePredictorRange, sep=""), collapse= " + ") #Creates Predictor String for Formula     
}

print("Status: Done")

# 
# Run SVM Models per Drug 
# 

print(paste("Status: Generating SVM Models & Predicting Using", nGenes, "Genes. . ."))
tPredictions        <- list()
for (i in 1:nDrugs) {
  knownClasses    <- vector()
  for (known in trainingData[drugs[i]]) knownClasses<- c(knownClasses, as.numeric(known))
  formula         <- as.formula(paste("factor(", drugs[i], " ~ ", predictorGenes, sep=""))
  model           <- svm(formula, trainingData, type = svmType, gamma = svmGamma, cost = svmCost, kernel = svmKernel, degree = svmDegree, coef0 = svmCoef0, cross = svmCross)
  predict         <- as.numeric(predict(model, trainingData))-1
  error           <- sum(trainingData[drugs[i]] - predict) / length(predict) * 100
  print(paste("SVM Model", i, "   Kernel:", svmKernel, "  Prediction Error:", decimals(abs(error), 2), "%     Drug:", drugs[i]))
  tPredicted      <- list(abs(as.numeric(predict(model, testingData)) - 1))
  tPredictions    <- cbind(tPredictions, tPredicted)
  AUC             <- prediction(predict, knownClasses)
  AUCPerf         <- performance(AUC, "tpr", "fpr")
  # plot(AUCPerf)
}
print("Status: Done")

#random Forests
library(randomForest)

bestmtry <- tuneRF(t(as.data.frame(trainExpressionData.mod.reorder)),as.factor(trainKey$Paclitaxel), ntreeTry=100, 
                   stepFactor=1.5,improve=0.01, trace=TRUE, plot=TRUE, dobest=FALSE)
rf1<-randomForest(formula, trainingData, ntree=999, importance=TRUE, keep.forest=TRUE)
j=2
nclass<-sum(trainKey[,j])
priorclass<-sum(trainKey[,j])/25
rf2<-randomForest(t(as.data.frame(trainExpressionData.mod.reorder)), 
                  y=as.factor(trainKey[,j]), mtry=91, ntree=999, keep.forest=TRUE, 
                  importance=TRUE, classwt=c(1-priorclass, priorclass))

adult.rf.pr = predict(adult.rf,type="prob",newdata=data$val)[,2]
adult.rf.pred = prediction(adult.rf.pr, data$val$income)
adult.rf.perf = performance(adult.rf.pred,"tpr","fpr")
plot(adult.rf.perf,main="ROC Curve for Random Forest",col=2,lwd=2)
abline(a=0,b=1,lwd=2,lty=2,col="gray")

importance(rf2)
varImpPlot(rf2)

MDSplot(rf2, trainKey[,j])
#hierarchical clustering
d=dist(as.matrix(t(expressionData[,-1])))
hc<-hclust(d)
plot(hc)
groups <- cutree(hc, k=4) # cut tree into 4 clusters
# draw dendogram with red borders around the 4 clusters 
rect.hclust(hc, k=4, border="red")

# Ward Hierarchical Clustering with Bootstrapped p values
library(pvclust)
fit <- pvclust(as.matrix(expressionData[,-1]), method.hclust="ward.d",
               method.dist="euclidean")
plot(fit) # dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit, alpha=.95)

# K-Means Clustering with 4 clusters
kmfit <- kmeans(t(expressionData[,-1]), 3)

# Cluster Plot against 1st 2 principal components

# vary parameters for most readable graph
library(cluster) 
#pammy<-pam(t(expressionData[,-1]), 4, metric="euclidean")
clusplot(t(expressionData[,-1]), kmfit$cluster, color=TRUE, shade=TRUE, labels=2, lines=0)

# Centroid Plot against 1st 2 discriminant functions
library(fpc)
plotcluster(t(expressionData[,-1]), kmfit$cluster)

expr.pca<-prcomp(t(expressionData[,-1]))
plot(expr.pca$x[,1], expr.pca$x[,2], col=(as.integer(kmfit$cluster)), pch=20, xlab="PCA1", ylab="PCA2", main="K-means clustering")
text(expr.pca$x[,1], expr.pca$x[,2], names(expressionData[,-1]), col=as.integer(kmfit$cluster))

#Test PCA
test.pca<-prcomp(t(testExpressionData.mod.reorder))
plot(test.pca$x[,1], test.pca$x[,2], pch=20, xlab="PCA1", ylab="PCA2", main="Test Set PCA")
text(test.pca$x[,1], test.pca$x[,2], names(testExpressionData.mod.reorder))