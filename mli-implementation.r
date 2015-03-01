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
# data(crabs)
# set.seed(1337)
# kp = sample(1:200, size=100)
# nn1 = MLearn(sp~CW+RW, data=crabs, nnetI, kp, size=100, decay=.020 )
# nn1

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
# Format Data for MLInterfaces - See Tutorial Here: http://www.bioconductor.org/packages//2.7/bioc/manuals/MLInterfaces/man/MLInterfaces.pdf
#

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
drugs               <- trainKeyTransposed$Drug # Vector of Drug Names
totalGeneCount      <- dim(as.data.frame(trainExpressionData))[1] # This is the Maximum Number of Predictors Possible
genePredictorRange  <- 1:100 #Selects Genes Labeled X1 - X10 in Training Data Set - 100 genes, Trying not to Get Too Crazy
predictorGenes      <- paste(paste("X", genePredictorRange, sep=""), collapse= " + ") #Creates Predictor String for Formula 

# MLInterfaces Common Parameters
set.seed(100) # random number generator, set seed specifies the seed so you can get consistent results
kp                  <- sample(1:100, size = 10) # samples from vector x (1:100) of size n (10)

#
# Test MLInterfaces Using Carboplatin Drug
#

# Generate Formula     
formula         <- as.formula(paste("Carboplatin ~ ", predictorGenes, sep=""))
# Neural Net
nn1             <- MLearn(formula, data = trainingData, nnetI, kp, size = 3, decay = .02, trace = TRUE)
# SVM
svm1            <- MLearn(formula, data = trainingData, svmI, kp)

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
