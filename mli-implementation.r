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

#We need to reorder the columns according to trainKey
trainExpressionData.mod.reorder<-trainExpressionData.mod[,testedCellLines]

# Mash Data Frames Together & Remove Duplicate Row Names From trainKey 
trainingData        <- data.frame(cbind(trainKey[-1], t(as.data.frame(trainExpressionData.mod.reorder))))
# And Testing Data - Naming Conventions (Not Sure if This Will Actually be Needed - Does MLInterfaces handle train/validate/test?)
testingData         <- testExpressionData

#
# Setup Parameters
#

# General Parameters
drugs               <- gsub("-", ".", trainKeyTransposed$Drug) # Vector of Drug Names
totalGeneCount      <- dim(as.data.frame(trainExpressionData))[1] # This is the Maximum Number of Predictors Possible
genePredictorRange  <- 1:10 #Selects Genes Labeled X1 - X10 in Training Data Set - 100 genes, Trying not to Get Too Crazy
predictorGenes      <- paste(paste("X", genePredictorRange, sep=""), collapse= " + ") #Creates Predictor String for Formula 

# MLInterfaces Common Parameters
set.seed(100) # random number generator, set seed specifies the seed so you can get consistent results
kp                  <- sample(1:25, size = 20) # samples from vector x (1:100) of size n (10)

#
# Test MLInterfaces Using Carboplatin Drug
#

# Generate Formula     
formula         <- as.formula(paste("Carboplatin ~ ", predictorGenes, sep=""))

# 
#  SVM Action
# 

# Run SVMs on Each Drug 

generateModels  <- function (response, predictors, df, t = "C-classification", g = 0.1, c = 100, k = "linear") {
    output <- list()
    for (r in response) {
        formula <- as.formula(paste(r, " ~ ", predictors, sep=""))
        model   <- svm(formula, data = df, type = t, gamma = g, cost = c, kernel = k)   
        output  <- c(output, model)
    }
    output
}

svmModels   <- generateModels(drugs, predictorGenes, trainingData)

# checkModels     <- function (response, models, df) {
#     for (model in models) {
#         prediction <- predict(model, )
#             p   <- predict(model, df)
#             trainLinearError        <- sum(trainingData$Carboplatin - p) / length(p) * 100
#             print(paste("Training Data SVM Linear Kernel Prediction Error:", decimals(abs(trainLinearError), 2), "%"))
#     }    
# }
# 

# 
# svmCheck    <- checkModels(drugs, svmModels, trainingData)

# checkModels <- function ()

# for (drug in drugs) {
#     # Generate Formula     
#     formula         <- as.formula(paste(drug, " ~ ", predictorGenes, sep=""))
#     # Neural Net
#     nn1             <- MLearn(formula, data = trainingData, nnetI, kp, size = 3, decay = .02, trace = TRUE)
#     # SVM
#     svm1            <- MLearn(formula, data = trainingData, svmI, kp)
# }

# # Generate SVM Models 
# linearModel             <- svm(formula, data = trainingData, Type = "c.classification", gamma = 00.1, cost = 100, kernel = "linear")
# print(linearModel)
# polyModel               <- svm(formula, data = trainingData, Type = "c.classification", gamma = 00.1, cost = 100, kernel = "polynomial")
# print(polyModel)
# radialModel             <- svm(formula, data = trainingData, Type = "c.classification", gamma = 00.1, cost = 100, kernel = "radial")
# print(radialModel)
# sigmoidModel            <- svm(formula, data = trainingData, Type = c.classification, gamma = 00.1, cost = 100, kernel = "sigmoid")
# print(sigmoidModel)

# # Linear Kernel Predictions 
# trainLinearPrediction   <- round(predict(linearModel, trainingData))
# trainLinearError        <- sum(trainingData$Carboplatin - trainLinearPrediction) / length(trainLinearPrediction) * 100
# print(paste("Training Data SVM Linear Kernel Prediction Error:", decimals(abs(trainLinearError), 2), "%"))

# testLinearPrediction    <- round(predict(linearModel, teData))
# testLinearError         <- sum(classTestData - testLinearPrediction) / length(testLinearPrediction) * 100
# print(paste("Testing Data SVM Linear Kernel Prediction Error:", decimals(abs(testLinearError), 2), "%"))

# # Polynomial Kernel Predictions 
# trainPolyPrediction     <- round(predict(polyModel, trainingData))
# trainPolyError          <- sum(trainingData$Carboplatin - trainPolyPrediction) / length(trainPolyPrediction) * 100
# print(paste("Training Data SVM Polynomial Kernel Prediction Error:", decimals(abs(trainPolyError), 2), "%"))

# testPolyPrediction      <- round(predict(polyModel, teData))
# testPolyError           <- sum(classTestData - testPolyPrediction) / length(testPolyPrediction) * 100
# print(paste("Testing Data SVM Polynomial Kernel Prediction Error:", decimals(abs(testPolyError), 2), "%"))

# # Radial Kernel Predictions 
# trainRadialPrediction   <- round(predict(radialModel, trainingData))
# trainRadialError        <- sum(trainingData$Carboplatin - trainRadialPrediction) / length(trainRadialPrediction) * 100
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