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

# setwd("~/Documents/Education/Graduate/OHSU/Courses/Winter 2015/Statistical Methods/assignments/project1/src") 

# 
# Utility Functions
# 

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

# General Parameters
drugs               <- gsub("-", ".", trainKeyTransposed$Drug) # Vector of Drug Names
totalGeneCount      <- dim(as.data.frame(trainExpressionData))[1] # This is the Maximum Number of Predictors Possible
genePredictorRange  <- 1:10 #Selects Genes Labeled X1 - X10 in Training Data Set - 100 genes, Trying not to Get Too Crazy
predictorGenes      <- paste(paste("X", genePredictorRange, sep=""), collapse= " + ") #Creates Predictor String for Formula 

# SVM Parameters 
svmCost     <- 100
svmGamma    <- 0.1
svmKernel   <- "linear"
svmType     <- "C-classification"

# 
# SVM Models per Drug 
# 

# SVM 1 - CGC.11047  
formula1    <- as.formula(paste("CGC.11047  ~ ", predictorGenes, sep=""))
model1      <- svm(formula1, trainingData, type = svmType, gamma = svmGamma, cost = svmCost, kernel = svmKernel)
predict1    <- predict(model1, trainingData)
error1      <- sum(trainingData$CGC.11047 - (as.numeric(predict1) - 1)) / length(predict1) * 100
print(paste("SVM Model 1 - Drug: CGC.11047 - Kernel:", svmKernel, "- Prediction Error:", decimals(abs(error1), 2), "%"))

# SVM 2 - Carboplatin    
formula2    <- as.formula(paste("Carboplatin ~ ", predictorGenes, sep=""))
model2      <- svm(formula2, trainingData, type = svmType, gamma = svmGamma, cost = svmCost, kernel = svmKernel)
predict2    <- predict(model2, trainingData)
error2      <- sum(trainingData$Carboplatin - (as.numeric(predict2) - 1)) / length(predict2) * 100
print(paste("SVM Model 2 - Drug: Carboplatin - Kernel:", svmKernel, "- Prediction Error:", decimals(abs(error2), 2), "%"))

# SVM 3 - Cisplatin    
formula3    <- as.formula(paste("Cisplatin ~ ", predictorGenes, sep=""))
model3      <- svm(formula3, trainingData, type = svmType, gamma = svmGamma, cost = svmCost, kernel = svmKernel)
predict3    <- predict(model3, trainingData)
error3      <- sum(trainingData$Cisplatin - (as.numeric(predict3) - 1)) / length(predict3) * 100
print(paste("SVM Model 3 - Drug: Cisplatin - Kernel:", svmKernel, "- Prediction Error:", decimals(abs(error3), 2), "%"))

# SVM 4 - GSK1070916  
formula4    <- as.formula(paste("GSK1070916 ~ ", predictorGenes, sep=""))
model4      <- svm(formula4, trainingData, type = svmType, gamma = svmGamma, cost = svmCost, kernel = svmKernel)
predict4    <- predict(model4, trainingData)
error4      <- sum(trainingData$GSK1070916 - (as.numeric(predict4) - 1)) / length(predict4) * 100
print(paste("SVM Model 4 - Drug: GSK1070916 - Kernel:", svmKernel, "- Prediction Error:", decimals(abs(error4), 2), "%"))

# SVM 5 - GSK1120212   
formula5    <- as.formula(paste("GSK1120212 ~ ", predictorGenes, sep=""))
model5      <- svm(formula5, trainingData, type = svmType, gamma = svmGamma, cost = svmCost, kernel = svmKernel)
predict5    <- predict(model5, trainingData)
error5      <- sum(trainingData$GSK1120212 - (as.numeric(predict5) - 1)) / length(predict5) * 100
print(paste("SVM Model 5 - Drug: GSK1120212 - Kernel:", svmKernel, "- Prediction Error:", decimals(abs(error5), 2), "%"))

# SVM 6 - GSK461364    
formula6    <- as.formula(paste("GSK461364 ~ ", predictorGenes, sep=""))
model6      <- svm(formula6, trainingData, type = svmType, gamma = svmGamma, cost = svmCost, kernel = svmKernel)
predict6    <- predict(model6, trainingData)
error6      <- sum(trainingData$GSK461364 - (as.numeric(predict6) - 1)) / length(predict6) * 100
print(paste("SVM Model 6 - Drug: GSK461364 - Kernel:", svmKernel, "- Prediction Error:", decimals(abs(error6), 2), "%"))

# SVM 7 - Geldanamycin   
formula7    <- as.formula(paste("Geldanamycin ~ ", predictorGenes, sep=""))
model7      <- svm(formula7, trainingData, type = svmType, gamma = svmGamma, cost = svmCost, kernel = svmKernel)
predict7    <- predict(model7, trainingData)
error7      <- sum(trainingData$Geldanamycin - (as.numeric(predict7) - 1)) / length(predict7) * 100
print(paste("SVM Model 7 - Drug: Geldanamycin - Kernel:", svmKernel, "- Prediction Error:", decimals(abs(error7), 2), "%"))

# SVM 8 - Oxaliplatin 
formula8    <- as.formula(paste("Oxaliplatin ~ ", predictorGenes, sep=""))
model8      <- svm(formula8, trainingData, type = svmType, gamma = svmGamma, cost = svmCost, kernel = svmKernel)
predict8    <- predict(model8, trainingData)
error8      <- sum(trainingData$Oxaliplatin - (as.numeric(predict8) - 1)) / length(predict8) * 100
print(paste("SVM Model 8 - Drug: Oxaliplatin - Kernel:", svmKernel, "- Prediction Error:", decimals(abs(error8), 2), "%"))

# SVM 9 - PF.3084014
formula9    <- as.formula(paste("PF.3084014 ~ ", predictorGenes, sep=""))
model9      <- svm(formula9, trainingData, type = svmType, gamma = svmGamma, cost = svmCost, kernel = svmKernel)
predict9    <- predict(model9, trainingData)
error9      <- sum(trainingData$PF.3084014 - (as.numeric(predict9) - 1)) / length(predict9) * 100
print(paste("SVM Model 9 - Drug: PF.3084014 - Kernel:", svmKernel, "- Prediction Error:", decimals(abs(error9), 2), "%"))

# SVM 10 - PF.3814735 
formula10    <- as.formula(paste("PF.3814735 ~ ", predictorGenes, sep=""))
model10      <- svm(formula10, trainingData, type = svmType, gamma = svmGamma, cost = svmCost, kernel = svmKernel)
predict10    <- predict(model10, trainingData)
error10      <- sum(trainingData$PF.3814735 - (as.numeric(predict10) - 1)) / length(predict10) * 100
print(paste("SVM Model 10 - Drug: PF.3814735 - Kernel:", svmKernel, "- Prediction Error:", decimals(abs(error10), 2), "%"))

# SVM 11 - PF.4691502 
formula11    <- as.formula(paste("PF.4691502 ~ ", predictorGenes, sep=""))
model11      <- svm(formula11, trainingData, type = svmType, gamma = svmGamma, cost = svmCost, kernel = svmKernel)
predict11    <- predict(model11, trainingData)
error11      <- sum(trainingData$PF.4691502 - (as.numeric(predict11) - 1)) / length(predict11) * 100
print(paste("SVM Model 11 - Drug: PF.4691502 - Kernel:", svmKernel, "- Prediction Error:", decimals(abs(error11), 2), "%"))

# SVM 12 - Paclitaxel
formula12    <- as.formula(paste("Paclitaxel ~ ", predictorGenes, sep=""))
model12      <- svm(formula12, trainingData, type = svmType, gamma = svmGamma, cost = svmCost, kernel = svmKernel)
predict12    <- predict(model12, trainingData)
error12      <- sum(trainingData$Paclitaxel - (as.numeric(predict12) - 1)) / length(predict12) * 100
print(paste("SVM Model 12 - Drug: Paclitaxel - Kernel:", svmKernel, "- Prediction Error:", decimals(abs(error12), 2), "%"))








############################################################################################################
# Old Code
############################################################################################################

# Run SVMs on Each Drug 

# generateModels  <- function (response, predictors, df, t = "C-classification", g = 0.1, c = 100, k = "linear") {
#     output <- list()
#     i = 1
#     for (r in response) {
#         formula <- as.formula(paste(r, " ~ ", predictors, sep=""))
#         model   <- svm(formula, data = df, type = t, gamma = g, cost = c, kernel = k)   
#         p       <- predict(model, df)
#         error   <- sum(df[i] - p) / length(p) * 100   
#         print(paste("Training Data SVM", k, "Kernel Prediction Error:", decimals(abs(error), 2), "%"))
#         output  <- c(output, model)
#         i = i + 1
#     }
#     output
# }

# svmModels   <- generateModels(drugs, predictorGenes, trainingData)

# 
# svmCheck    <- checkModels(drugs, svmModels, trainingData)

# checkModels <- function ()

# # Generate SVM Models 
# linearModel             <- svm(formula, data = trainingData, type = "C-classification", gamma = 00.1, cost = 100, kernel = "linear")
# print(linearModel)
# polyModel               <- svm(formula, data = trainingData, type = "C-classification", gamma = 00.1, cost = 100, kernel = "polynomial")
# print(polyModel)
# radialModel             <- svm(formula, data = trainingData, type = "C-classification", gamma = 00.1, cost = 100, kernel = "radial")
# print(radialModel)
# sigmoidModel            <- svm(formula, data = trainingData, type = "C-classification", gamma = 00.1, cost = 100, kernel = "sigmoid")
# print(sigmoidModel)

# Linear Kernel Predictions 
#trainLinearPrediction   <- round(predict(linearModel, trainingData))
# trainLinearPrediction   <- predict(linearModel, trainingData)
# trainLinearError        <- sum(trainingData$Carboplatin - (as.numeric(trainLinearPrediction)-1)) / length(trainLinearPrediction) * 100
# print(paste("Training Data SVM Linear Kernel Prediction Error:", decimals(abs(trainLinearError), 2), "%"))
# 
# # testLinearPrediction    <- round(predict(linearModel, teData))
# # testLinearError         <- sum(classTestData - testLinearPrediction) / length(testLinearPrediction) * 100
# # print(paste("Testing Data SVM Linear Kernel Prediction Error:", decimals(abs(testLinearError), 2), "%"))
# 
# # Polynomial Kernel Predictions 
# #trainPolyPrediction     <- round(predict(polyModel, trainingData))
# trainPolyPrediction     <- predict(polyModel, trainingData)
# trainPolyError          <- sum(trainingData$Carboplatin - (as.numeric(trainPolyPrediction)-1)) / length(trainPolyPrediction) * 100
# print(paste("Training Data SVM Polynomial Kernel Prediction Error:", decimals(abs(trainPolyError), 2), "%"))
# 
# # testPolyPrediction      <- round(predict(polyModel, teData))
# # testPolyError           <- sum(classTestData - testPolyPrediction) / length(testPolyPrediction) * 100
# # print(paste("Testing Data SVM Polynomial Kernel Prediction Error:", decimals(abs(testPolyError), 2), "%"))
# 
# # Radial Kernel Predictions 
# #trainRadialPrediction   <- round(predict(radialModel, trainingData))
# trainRadialPrediction   <- predict(radialModel, trainingData)
# trainRadialError        <- sum(trainingData$Carboplatin - (as.numeric(trainRadialPrediction)-1)) / length(trainRadialPrediction) * 100
# print(paste("Training Data SVM Radial Kernel Prediction Error:", decimals(abs(trainRadialError), 2), "%"))

# testRadialPrediction    <- round(predict(radialModel, teData))
# testRadialError         <- sum(classTestData - testRadialPrediction) / length(testRadialPrediction) * 100
# print(paste("Testing Data SVM Radial Kernel Prediction Error:", decimals(abs(testRadialError), 2), "%"))

# Sigmoid Kernel Predictions 
# trainSigmoidPrediction  <- round(predict(sigmoidModel, trainingData))
# trainSigmoidError       <- sum(trainingData$Carboplatin - trainSigmoidPrediction) / length(trainSigmoidPrediction) * 100
# print(paste("Training Data SVM Sigmoid Kernel Prediction Error:", decimals(abs(trainSigmoidError), 2), "%"))

# testSigmoidPrediction   <- round(predict(sigmoidModel, teData))
# testSigmoidError        <- sum(classTestData - testSigmoidPrediction) / length(testSigmoidPrediction) * 100
# print(paste("Testing Data SVM Sigmoid Kernel Prediction Error:", decimals(abs(testSigmoidError), 2), "%"))

#
# ML Interfaces Stuff
#

# Neural Net
# nn1             <- MLearn(formula, data = trainingData, nnetI, kp, size = 10, decay = .02, trace = TRUE)
# nn1
# SVM
# svm1            <- MLearn(formula, data = trainingData, svmI, kp)
# svm1

#
# Run MLInterfaces Routine on All Drugs
#

# Run ML Algorithms on Each Drug 
# for (drug in drugs) {
#     # Generate Formula     
#     formula         <- as.formula(paste(drug, " ~ ", predictorGenes, sep=""))
#     # Neural Net
#     nn1             <- MLearn(formula, data = trainingData, nnetI, kp, size = 3, decay = .02, trace = TRUE)
#     # SVM
#     svm1            <- MLearn(formula, data = trainingData, svmI, kp)
# }