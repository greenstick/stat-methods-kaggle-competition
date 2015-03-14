#!/usr/bin/env Rscript

#
# Setup
#

print("Status: Setup")

# Installations
source("http://bioconductor.org/biocLite.R")
# biocLite("MLInterfaces")
# install.packages("MLInterfaces")
# install.packages('randomForest')
# install.packages("ROCR")
# library("MASS")
library(ROCR)
library(e1071)

# Get Session Information
# sessionInfo()

# Set Directory
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
    y <- format(round(x, d), nsmall = d)
    y
}

# Get Coefficient of Variation
cv                  <- function (x, dimension = 2) {
    if (is.data.frame(x) == TRUE) {
        y <- 100 * (apply(x, dimension, sd, na.rm=TRUE) / apply(x, dimension, mean, na.rm=TRUE))
    } else {
        y <- 100 * (sd(x, na.rm=TRUE) / mean(x, na.rm=TRUE)) 
    }
    y
}

#Not In function
`%notin%` <- function (x,y) !(x %in% y) 

# Gathers high CV genes and targets and reassigns them to training data
getByCV         <- function (df, vThreshold = 2, geneStartColumn = 1, absoluteVal = TRUE, overThresh = TRUE) {
    v           <- vector()
    for (i in geneStartColumn:ncol(df)) {
        l               <- length(v)
        cVar            <- cv(df[i])
        v               <- c(v, cVar)
    }
    if (absoluteVal == TRUE) {
        v <- which(abs(t(as.matrix(v))) > vThreshold)
    } else {
        if (overThresh == TRUE) {
            v <- which(t(as.matrix(v)) > vThreshold)
        } else {
            v <- which(t(as.matrix(v)) < vThreshold)
        }
    }
    v 
}

getByMax            <- function (df, mThreshold = 4, geneStartColumn = 1, absoluteVal) {
    v <- vector()
    for (i in 1:ncol(df)) {
        l <- length(v)
        if (absoluteVal == TRUE) {
            cMax <- as.numeric(max(abs(df[i])))
        } else {
            cMax <- as.numeric(max(df[i]))
        }
        v <- x(v, cMax)
    }
    v <- which(t(as.matrix(v)) > mThreshold)
    v
}

print("Status: Done")

#
# Data Imports
# 

print("Status: Importing Data . . .")

subtypeData         <- read.csv("data/subtype.csv")
mapData             <- read.csv("data/scoring_and_test_set_id_mappings.csv")
trainKey            <- read.csv("data/training_set_answers.csv")
trainKeyTransposed  <- read.csv("data/training_set_answers_transposed.csv")
expressionData      <- read.csv("data/expression.csv")

print("Status: Done")

#
# Data Wrangling
# 

print("Status: Configuring . . .")

# Format Data for MLInterfaces - See Tutorial Here: http://www.bioconductor.org/packages//2.7/bioc/manuals/MLInterfaces/man/MLInterfaces.pdf
testedCellLines     <- (trainKey[1])$X

predictCellLines    <- c("HCC1187", "MCF7", "MDAMB361", "MDAMB231", "BT549", "X600MPE", "HCC1954", "SKBR3", "MCF10A", "MCF12A", "HCC3153", "MDAMB157", "LY2", "AU565")
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
testedCellLines.mod[4]<-paste("X", testedCellLines.mod[4], sep="") #changes 184A1 to X184A1

#Extracting desired cell line columns out of expressionData
trainExpressionData.mod<-expressionData[which(names(expressionData) %in% testedCellLines.mod)]
testExpressionData.mod<-expressionData[which(names(expressionData) %notin% testedCellLines.mod)]

#We need to reorder the columns according to trainKey
trainExpressionData.mod.reorder<-trainExpressionData.mod[, testedCellLines]
testExpressionData.mod.reorder<-testExpressionData.mod[, predictCellLines]

# Mash Data Frames Together & Remove Duplicate Row Names From trainKey 
trainingData        <- data.frame(cbind(trainKey[-1], t(as.data.frame(trainExpressionData.mod.reorder))))
# And Testing Data - Naming Conventions (Not Sure if This Will Actually be Needed - Does MLInterfaces handle train/validate/test?)
testingData         <- data.frame(t(as.data.frame(testExpressionData.mod.reorder)))

print("Status: Done")

#
# Setup Parameters
#

print("Status: Loading parameters . . .")

# SVM Parameters 
svmCost     <- 20
svmGamma    <- 0.078
svmKernel   <- "linear"
svmDegree   <- 2
svmType     <- "one-classification"
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
    print(paste("Generating SVM Model", i, "& Predicting . . ."))
    knownClasses    <- vector()
    for (known in trainingData[drugs[i]]) knownClasses<- c(knownClasses, as.numeric(known))
    formula         <- as.formula(paste(drugs[i], "~", predictorGenes))
    model           <- svm(formula, trainingData, type = svmType, gamma = svmGamma, cost = svmCost, kernel = svmKernel, degree = svmDegree, coef0 = svmCoef0, cross = svmCross)
    predict         <- as.numeric(predict(model, trainingData)) - 1
    error           <- sum(trainingData[drugs[i]] - predict) / length(predict) * 100
    print(paste("SVM Model", i, "   Kernel:", svmKernel, "  Prediction Error:", decimals(abs(error), 2), "%     Drug:", drugs[i]))
    tPredicted      <- list(abs(as.numeric(predict(model, testingData)) - 1))
    tPredictions    <- cbind(tPredictions, tPredicted)
    AUC             <- prediction(predict, knownClasses)
    AUCPerf         <- performance(AUC, "tpr")
    # plot(AUCPerf)
}
print("Status: Done")

#
# Wrangle Output & Save Predictions to CSV
#

print("Status: Generating Output File")
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
filename     <- paste("submissions/svm/genes_", nGenes, "_cvSelect_", SelectCV , "_kernel_", svmKernel, "_cost_", svmCost, "_gamma_", svmGamma, "_degree_", svmDegree, "_coef0_", svmCoef0, "_cross_", svmCross, ".csv", sep="")
write.csv(submission, file=filename, row.names = FALSE)
print(paste("Output Saved As:", filename))

print("Status: Done")