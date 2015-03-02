#!/usr/bin/env Rscript

#
# Setup
#

# Installations
# source("http://bioconductor.org/biocLite.R")
# biocLite("MLInterfaces")
# install.packages("MLInterfaces")
# library("MASS")
library(MLInterfaces)
library(e1071)

setwd("~/Documents/Education/Graduate/OHSU/Courses/Winter 2015/Statistical Methods/assignments/project1/src") 

# 
# Utility Functions
# 

# Insert Into Vector 
insert              <- function (vector, element, position = 0) {
    length <- length(vector)
    if (position == 0) {
        return (c(element, vector[0:(length)]))
    } else if (position > length) {
        return (c(vector[1:(length)], element))       
    } else {
        return (c(vector[1:(position - 1)], element, vector[(position):length]))
    }
}

# How Many Decimals?
decimals            <- function (x, d) {
    x <- format(round(x, d), nsmall = d)
    x
}

#Not In function
`%notin%` <- function (x,y) !(x %in% y) 

#
# Data Imports
# 

# setwd Set Working Directory to Kaggle repository

subtypeData         <- read.csv("data/subtype.csv")
mapData             <- read.csv("data/scoring_and_test_set_id_mappings.csv")
trainKey            <- read.csv("data/training_set_answers.csv")
trainKeyTransposed  <- read.csv("data/training_set_answers_transposed.csv")
expressionData      <- read.csv("data/expression.csv")

#
# Data Wrangling
# 

# Format Data for MLInterfaces - See Tutorial Here: http://www.bioconductor.org/packages//2.7/bioc/manuals/MLInterfaces/man/MLInterfaces.pdf
testedCellLines     <- (trainKey[1])$X
predictCellLines <- c("HCC1187", "MCF7", "MDAMB361", "MDAMB231", "BT549", "X600MPE", "HCC1954", "SKBR3", "MCF10A", "MCF12A", "HCC3153", "MDAMB157", "LY2", "AU565")
trainExpressionData <- list()
testExpressionData  <- list()
for (value in testedCellLines) {
    for (j in 2:ncol(expressionData)) {
        if (value == names(expressionData[j]) || paste("X", value, sep="") == names(expressionData[j])) {
            trainExpressionData <- c(trainExpressionData, expressionData[j])
        } else {
            testExpressionData  <- c(testExpressionData, expressionData[j])  
        }
    }
}

testedCellLines.mod<-lapply(trainKey$X, as.character)
testedCellLines.mod[4]<-paste("X", testedCellLines.mod[4],sep="") #changes 184A1 to X184A1

#Extracting desired cell line columns out of expressionData
trainExpressionData.mod<-expressionData[which(names(expressionData) %in% testedCellLines.mod)]
testExpressionData.mod<-expressionData[which(names(expressionData) %notin% testedCellLines.mod)]

#We need to reorder the columns according to trainKey
trainExpressionData.mod.reorder<-trainExpressionData.mod[,testedCellLines]
testExpressionData.mod.reorder<-testExpressionData.mod[, predictCellLines]

# Mash Data Frames Together & Remove Duplicate Row Names From trainKey 
trainingData        <- data.frame(cbind(trainKey[-1], t(as.data.frame(trainExpressionData.mod.reorder))))
# And Testing Data - Naming Conventions (Not Sure if This Will Actually be Needed - Does MLInterfaces handle train/validate/test?)
testingData         <- data.frame(t(as.data.frame(testExpressionData.mod.reorder)))

#
# Setup Parameters
#

# SVM Parameters 
svmCost     <- 1
svmGamma    <- 0.12
svmKernel   <- "polynomial"
svmDegree   <- 3
svmType     <- "C-classification"
svmCoef0    <- 2
svmCross    <- 1
nStart      <- 1
nGenes      <- 500

# General Parameters
drugs               <- gsub("-", ".", trainKeyTransposed$Drug) # Vector of Drug Names
totalGeneCount      <- dim(as.data.frame(trainExpressionData))[1] # This is the Maximum Number of Predictors Possible
genePredictorRange  <- nStart:(nStart + nGenes - 1) #Selects Genes Labeled X1 - X10 in Training Data Set - 100 genes, Trying not to Get Too Crazy
predictorGenes      <- paste(paste("X", genePredictorRange, sep=""), collapse= " + ") #Creates Predictor String for Formula 
print(genePredictorRange[2])

# 
# SVM Models per Drug 
# 

# SVM 1 - CGC.11047  
formula1    <- as.formula(paste("CGC.11047  ~ ", predictorGenes, sep=""))
model1      <- svm(formula1, trainingData, type = svmType, gamma = svmGamma, cost = svmCost, kernel = svmKernel, degree = svmDegree, coef0 = svmCoef0, cross = svmCross)
predict1    <- predict(model1, trainingData)
error1      <- sum(trainingData$CGC.11047 - (as.numeric(predict1) - 1)) / length(predict1) * 100
print(paste("SVM Model 1 - Drug: CGC.11047 - Kernel:", svmKernel, "- Prediction Error:", decimals(abs(error1), 2), "%"))
tPredict1   <- as.numeric(predict(model1, testingData)) - 1

# SVM 2 - Carboplatin    
formula2    <- as.formula(paste("Carboplatin ~ ", predictorGenes, sep=""))
model2      <- svm(formula2, trainingData, type = svmType, gamma = svmGamma, cost = svmCost, kernel = svmKernel, degree = svmDegree, coef0 = svmCoef0, cross = svmCross)
predict2    <- predict(model2, trainingData)
error2      <- sum(trainingData$Carboplatin - (as.numeric(predict2) - 1)) / length(predict2) * 100
print(paste("SVM Model 2 - Drug: Carboplatin - Kernel:", svmKernel, "- Prediction Error:", decimals(abs(error2), 2), "%"))
tPredict2   <- as.numeric(predict(model2, testingData)) -1

# SVM 3 - Cisplatin    
formula3    <- as.formula(paste("Cisplatin ~ ", predictorGenes, sep=""))
model3      <- svm(formula3, trainingData, type = svmType, gamma = svmGamma, cost = svmCost, kernel = svmKernel, degree = svmDegree, coef0 = svmCoef0, cross = svmCross)
predict3    <- predict(model3, trainingData)
error3      <- sum(trainingData$Cisplatin - (as.numeric(predict3) - 1)) / length(predict3) * 100
print(paste("SVM Model 3 - Drug: Cisplatin - Kernel:", svmKernel, "- Prediction Error:", decimals(abs(error3), 2), "%"))
tPredict3   <- as.numeric(predict(model3, testingData)) - 1

# SVM 4 - GSK1070916  
formula4    <- as.formula(paste("GSK1070916 ~ ", predictorGenes, sep=""))
model4      <- svm(formula4, trainingData, type = svmType, gamma = svmGamma, cost = svmCost, kernel = svmKernel, degree = svmDegree, coef0 = svmCoef0, cross = svmCross)
predict4    <- predict(model4, trainingData)
error4      <- sum(trainingData$GSK1070916 - (as.numeric(predict4) - 1)) / length(predict4) * 100
print(paste("SVM Model 4 - Drug: GSK1070916 - Kernel:", svmKernel, "- Prediction Error:", decimals(abs(error4), 2), "%"))
tPredict4   <- as.numeric(predict(model4, testingData)) - 1

# SVM 5 - GSK1120212   
formula5    <- as.formula(paste("GSK1120212 ~ ", predictorGenes, sep=""))
model5      <- svm(formula5, trainingData, type = svmType, gamma = svmGamma, cost = svmCost, kernel = svmKernel, degree = svmDegree, coef0 = svmCoef0, cross = svmCross)
predict5    <- predict(model5, trainingData)
error5      <- sum(trainingData$GSK1120212 - (as.numeric(predict5) - 1)) / length(predict5) * 100
print(paste("SVM Model 5 - Drug: GSK1120212 - Kernel:", svmKernel, "- Prediction Error:", decimals(abs(error5), 2), "%"))
tPredict5   <- as.numeric(predict(model5, testingData)) - 1

# SVM 6 - GSK461364    
formula6    <- as.formula(paste("GSK461364 ~ ", predictorGenes, sep=""))
model6      <- svm(formula6, trainingData, type = svmType, gamma = svmGamma, cost = svmCost, kernel = svmKernel, degree = svmDegree, coef0 = svmCoef0, cross = svmCross)
predict6    <- predict(model6, trainingData)
error6      <- sum(trainingData$GSK461364 - (as.numeric(predict6) - 1)) / length(predict6) * 100
print(paste("SVM Model 6 - Drug: GSK461364 - Kernel:", svmKernel, "- Prediction Error:", decimals(abs(error6), 2), "%"))
tPredict6   <- as.numeric(predict(model6, testingData)) - 1

# SVM 7 - Geldanamycin   
formula7    <- as.formula(paste("Geldanamycin ~ ", predictorGenes, sep=""))
model7      <- svm(formula7, trainingData, type = svmType, gamma = svmGamma, cost = svmCost, kernel = svmKernel, degree = svmDegree, coef0 = svmCoef0, cross = svmCross)
predict7    <- predict(model7, trainingData)
error7      <- sum(trainingData$Geldanamycin - (as.numeric(predict7) - 1)) / length(predict7) * 100
print(paste("SVM Model 7 - Drug: Geldanamycin - Kernel:", svmKernel, "- Prediction Error:", decimals(abs(error7), 2), "%"))
tPredict7   <- as.numeric(predict(model7, testingData)) - 1

# SVM 8 - Oxaliplatin 
formula8    <- as.formula(paste("Oxaliplatin ~ ", predictorGenes, sep=""))
model8      <- svm(formula8, trainingData, type = svmType, gamma = svmGamma, cost = svmCost, kernel = svmKernel, degree = svmDegree, coef0 = svmCoef0, cross = svmCross)
predict8    <- predict(model8, trainingData)
error8      <- sum(trainingData$Oxaliplatin - (as.numeric(predict8) - 1)) / length(predict8) * 100
print(paste("SVM Model 8 - Drug: Oxaliplatin - Kernel:", svmKernel, "- Prediction Error:", decimals(abs(error8), 2), "%"))
tPredict8   <- as.numeric(predict(model8, testingData)) - 1

# SVM 9 - PF.3084014
formula9    <- as.formula(paste("PF.3084014 ~ ", predictorGenes, sep=""))
model9      <- svm(formula9, trainingData, type = svmType, gamma = svmGamma, cost = svmCost, kernel = svmKernel, degree = svmDegree, coef0 = svmCoef0, cross = svmCross)
predict9    <- predict(model9, trainingData)
error9      <- sum(trainingData$PF.3084014 - (as.numeric(predict9) - 1)) / length(predict9) * 100
print(paste("SVM Model 9 - Drug: PF.3084014 - Kernel:", svmKernel, "- Prediction Error:", decimals(abs(error9), 2), "%"))
tPredict9   <- as.numeric(predict(model9, testingData)) - 1

# SVM 10 - PF.3814735 
formula10    <- as.formula(paste("PF.3814735 ~ ", predictorGenes, sep=""))
model10      <- svm(formula10, trainingData, type = svmType, gamma = svmGamma, cost = svmCost, kernel = svmKernel, degree = svmDegree, coef0 = svmCoef0, cross = svmCross)
predict10    <- predict(model10, trainingData)
error10      <- sum(trainingData$PF.3814735 - (as.numeric(predict10) - 1)) / length(predict10) * 100
print(paste("SVM Model 10 - Drug: PF.3814735 - Kernel:", svmKernel, "- Prediction Error:", decimals(abs(error10), 2), "%"))
tPredict10   <- as.numeric(predict(model10, testingData)) - 1

# SVM 11 - PF.4691502 
formula11    <- as.formula(paste("PF.4691502 ~ ", predictorGenes, sep=""))
model11      <- svm(formula11, trainingData, type = svmType, gamma = svmGamma, cost = svmCost, kernel = svmKernel, degree = svmDegree, coef0 = svmCoef0, cross = svmCross)
predict11    <- predict(model11, trainingData)
error11      <- sum(trainingData$PF.4691502 - (as.numeric(predict11) - 1)) / length(predict11) * 100
print(paste("SVM Model 11 - Drug: PF.4691502 - Kernel:", svmKernel, "- Prediction Error:", decimals(abs(error11), 2), "%"))
tPredict11   <- as.numeric(predict(model11, testingData)) - 1

# SVM 12 - Paclitaxel
formula12    <- as.formula(paste("Paclitaxel ~ ", predictorGenes, sep=""))
model12      <- svm(formula12, trainingData, type = svmType, gamma = svmGamma, cost = svmCost, kernel = svmKernel, degree = svmDegree, coef0 = svmCoef0, cross = svmCross)
predict12    <- predict(model12, trainingData)
error12      <- sum(trainingData$Paclitaxel - (as.numeric(predict12) - 1)) / length(predict12) * 100
print(paste("SVM Model 12 - Drug: Paclitaxel - Kernel:", svmKernel, "- Prediction Error:", decimals(abs(error12), 2), "%"))
tPredict12   <- as.numeric(predict(model12, testingData)) - 1

#
# Wrangle Output & Save Predictions to CSV
#

tPredictions <- list(tPredict1, tPredict2, tPredict3, tPredict4, tPredict5, tPredict6, tPredict7, tPredict8, tPredict9, tPredict10, tPredict11, tPredict12)
subTop       <- list()
subBottom    <- list()
for (col in tPredictions) {
    top         <- col[1:9]
    subTop      <- insert(subTop, top, length(subTop))
    bottom      <- col[10:14]
    subBottom   <- insert(subBottom, bottom, length(subBottom))
}
values       <- t(cbind(t(subTop), t(subBottom)))
ids          <- seq(1, length(values), 1)
subdf        <- data.frame(id=as.matrix(ids), value=as.matrix(values))
submission   <- data.frame(lapply(subdf, as.character), stringsAsFactors=FALSE)
print(submission)
filename     <- paste("submissions/svm/genes_", nGenes, "_kernel_", svmKernel, "_cost_", svmCost, "_gamma_", svmGamma, "_degree_", svmDegree, "_coef0_", svmCoef0, "_cross_", svmCross, ".csv", sep="")
write.csv(submission, file=filename, row.names = FALSE)

