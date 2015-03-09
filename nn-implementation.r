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
library(ROCR)

# Get Session Information
# sessionInfo()

# Set Directory
setwd("~/Documents/Education/Graduate/OHSU/Courses/Winter 2015/Statistical Methods/assignments/project1/src") 

# Set Seed
set.seed(1337)

# 
# Some Utility Functions
# 

# Utility Function to Insert Element into Vector at Specified Position (position 0 default) 
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

# Utility Function to Retrieve Column From Data Frame - Just Because The Standard Syntax Gives Me a Headache
getColumn             <- function (df, col) {
    c <- (df[[col]])
    c
}

# Subsets a Data Frame Column by It's Name Using a Regex
getColumnByRegex  <- function (df, regex) {
    d <- subset(df, select=(names(df)[grep(regex, names(df))]))
    d
}

# Get Test Data - Cause Names And Data Are in Different Columns. LOL! 
getTestData         <- function (data, classNames, class) {
    index           <- 1
    l               <- list()
    for (name in classNames) {
        if (name == class) { 
            l <- c(l, data[index]) 
        }
        index       <- index + 1
    }
    l
}

# Coefficient of Variation
cv                  <- function (x) {
    y <- 100 * (sd(x, na.rm=TRUE) / mean(x, na.rm=TRUE)) 
    y
}

# Subsets matrix by columns by quantiles
getByQuantile   <- function (m, lowerCutoff = 0.25, upperCutoff = 0.75) {
    q           <- quantile(m, probs = c(lowerCutoff, upperCutoff))
    columns     <- apply(m, 2, function (x) any(x < q[1] | x > q[2]))
    m           <- m[ ,columns]
    m
}

# Subsets matrix by variance cutoff
getByCV   <- function (m, minCV) {
    columns     <- apply(m, 2, function (x) any(cv(x) > minCV))
    m           <- m[ ,columns]
    output      <- list(matrix = m, columns = columns)
    output
}

# Generates a Numeric Vector of Length n With a Specified Value or Random Values (default behavior)
generateVector      <- function (n, value = FALSE, lower = -0.5, upper = 0.5) {
    if (is.numeric(value == FALSE)) {
        print("Error: generateVector() - Value Should Be a Numeric Type")
        v = FALSE
    } else if (identical(value, FALSE)) {
        v <- c(runif(n, min = lower, max = upper))
    } else {
        v <- c(rep(value, n))
    }
    v
}

# How Mayn Decimals?
decimals            <- function (x, d) {
    y <- format(round(x, d), nsmall = d)
    y
}

# 
# Learning Kernels
# 

# Signum function -- Ensures Proper Evaluation of 0 Input 
signum              <- function (x) {
    s <- sign(x)
    s[s == 0] <- 1
    s
}

# Sigmoid function -- To Learn Good 
sigmoid             <- function (x) {
    y <- (1 / (1 + exp(-x)))
    y
}

# Fast Sigmoid Approximation
fSigmoid            <- function (x) {
    y <- x / (1 + abs(x))
    y
}

# 
# Artifical Neural Network Training Function
# 

ANN.train <- function (epochs = 2, inputTrainingData, trainingTargets, etaP = 0.1, etaH = 0.01, hiddenNodes = 20, dataWeightsLimit = 0.05, hiddenWeightsLimit = 0.5, plotData = list(), visualizeWeights = FALSE, verboseTrain = TRUE) {
    # Configuration     
    inputTrainingDataDim    <- dim(inputTrainingData)
    inputTrainingLengthC    <- inputTrainingDataDim[1]
    inputTrainingLengthR    <- inputTrainingDataDim[2]
    trainingWeights         <- matrix(runif(((inputTrainingLengthR) * hiddenNodes), min = - dataWeightsLimit, max = dataWeightsLimit), nrow = inputTrainingLengthR, ncol = hiddenNodes, byrow = T)
    hiddenWeights           <- matrix(generateVector(hiddenNodes + 1, lower = - hiddenWeightsLimit, upper = hiddenWeightsLimit), nrow = hiddenNodes + 1, ncol = 1, byrow = T)
    geneExpression          <- data.frame()
    if (length(plotData) > 0) {
        if (plotData$SSE == TRUE) {
            par(mfrow = c(1, 1))
            plot(x = NULL, y = NULL, xlim = c(1, epochs), ylim = c(0, 1), cex.main = 0.8, ylab = "Classification SSE", xlab = "Epoch", main = "Scatter Plot of Sample Classification \n SSE by Epoch")
        } else if (plotData$distance == TRUE) {
            par(mfrow = c(1, 1))
            plot(x = NULL, y = NULL, xlim = c(1, epochs), ylim = c(0, 1), cex.main = 0.8, ylab = "Classification Error Distance", xlab="Epoch", main = "Scatter Plot of Sample Classification \n Error Distance by Epoch")
        }
    }
    # Epoch Loop   
    if (verboseTrain == FALSE) {
        print("=========================================")
        print("-- Training . . .    --------------------")
        print("=========================================")  
    }
    for (n in 1:epochs) {
        # Configure Epoch
        distanceList        <- vector()
        iRandomization      <- sample(inputTrainingLengthC)
        hit                 <- 0
        miss                <- 0
        grand               <- 0
        meanSSE             <- vector()
        meanDistance        <- vector()
        classifications     <- vector()
        if (verboseTrain == TRUE) {
            print("=========================================")
            print(paste("-- Epoch", n, "Start"))
            print("=========================================")
        }
        # Sample Loop         
        for (i in iRandomization) {
            grand           <- grand + 1
            total           <- i * n
            if (verboseTrain == TRUE) {
                print("-- New Sample ---------------------------")
                print("   Forward-Propagating...")
                print(paste("       Sample", i, "of Epoch", n))
                print(paste("       Known Class:           ", trainingTargets[i]))
            }
            # Forward Propagation
            hiddenLayer     <- sigmoid(inputTrainingData[i,] %*% trainingWeights)
            hiddenLayer     <- c(insert(hiddenLayer, element = 1)) # Insert Bias Term to Hidden Layer
            output          <- sigmoid(hiddenLayer %*% hiddenWeights)
            classifications <- c(classifications, output)
            # Metric Computation & Communication   
            error           <- 0.5 * (trainingTargets[i] - output) ^ 2
            distance        <- abs(trainingTargets[i] - output) # An Easily Interpretable Error Measure
            distanceList    <- c(distanceList, distance) # For Monitoring Error Reduction
            rounded         <- round(output)
            if (verboseTrain == TRUE) {
                print(paste("       Computed Class:        ", decimals(output, 8)))
                print(paste("       Raw Distance:          ", decimals(distance, 8)))
                print(paste("       Computed SSE:          ", decimals(error, 8)))
            }
            if (abs(rounded - trainingTargets[i]) == 0) {
                hit <- hit + 1
                if (verboseTrain == TRUE) print("       Rounded Hit Status:     hit! :)")
            } else {
                miss <- miss + 1
                if (verboseTrain == TRUE) print("       Rounded Hit Status:     miss :(")
            }
            if (verboseTrain == TRUE) {
                print(paste("           Epoch Hits / Total:", hit, "/", grand))
                print(paste("           Epoch Hit Percent: ", decimals((hit/grand) * 100, 2)))
                print("   Back-Propagating...")
            }
            # Back Propagation
            deltaP          <- drop((trainingTargets[i] - output) * output * (1 - output))
            deltaH          <- (hiddenLayer * (1 - hiddenLayer) * hiddenWeights * deltaP)[-c(1)] # Compute deltaH & Remove Bias Term
            hiddenChange    <- etaH * hiddenLayer * deltaP
            updatedHidden   <- as.matrix(hiddenLayer) + as.matrix(hiddenChange)
            hiddenLayer     <- updatedHidden
            weightChange    <- etaP * t(inputTrainingData[i,, drop = FALSE]) %*% deltaH
            updatedWeights  <- trainingWeights + weightChange
            trainingWeights <- updatedWeights
            # Plot Points & Communication
            if (length(plotData) > 0) {
                if (plotData$SSE == TRUE) {
                    meanSSE     <- c(insert(meanSSE, element = error)) 
                    points(x = n, y = error , col = 'orange')
                } else if (plotData$distance == TRUE) {
                    meanDistance<- c(insert(meanDistance, element = distance))
                    points(x = n, y = distance, col = 'blue')   
                }
            }
            if (verboseTrain == TRUE) print("-- Sample Done --------------------------")
        }
        if (verboseTrain == TRUE) {
            print("=========================================")
            print(paste("-- Epoch", n, "Done"))
            print("=========================================")
        }
        if (length(plotData) > 0) {
            if (plotData$SSE == TRUE) {
                points(x = n, y = mean(meanSSE) , col = 'green')
            } else if (plotData$distance) {
                points(x = n, y = mean(meanDistance), col = 'red')   
            }
        }
    }
    if (verboseTrain == FALSE) {
        print("=========================================")
        print("-- Training Complete     ----------------")
        print("=========================================") 
        print("-- Report -------------------------------")
        print("   Rounded Hits Last Epoch:")
        print(paste("       Train Hits / Total:     ", hit, "/", grand))
        print(paste("       Train Hit Percent:      ", decimals((hit/grand) * 100, 2))) 
    }
    if (length(plotData) > 0) {
        if (plotData$SSE == TRUE) {
            legend(1, 1, legend = c('Mean SSE', 'SSE per Sample'), pch = 1, col = c("green", "orange"), cex = 0.6)   
        } else if (plotData$distance == TRUE) {
            legend(1, 1, legend = c('Mean Error Distance', 'Error Distance per Sample'), pch = 1, col = c("red", "blue"), cex = 0.6)   
        } else if (plotData$weightMeans$plot == TRUE) {
            geneExpression <- as.data.frame(t(apply(trainingWeights, 1, mean)))[-c(1)]
            names(geneExpression) <- geneNames
            image(t(as.matrix(geneExpression)), axis = FALSE, main = "Heat Map of Mean Weights Computed for \n Gene Expression Levels", axes = FALSE)
        }
    }
    # Return Results
    list (
        classifications     = classifications,
        trainedWeights      = trainingWeights, 
        trainedHidden       = hiddenWeights, 
        meanGeneExpression  = geneExpression
    )
}

# 
# Artifical Neural Network Classification Function
# 

ANN.classify <- function (inputTestData, testTargets = vector(), calibratedWeights = computedWeights$trainedWeights, calibratedHidden = computedWeights$trainedHidden, verboseClassify = TRUE) {
    # Configure
    inputTestDataDim    <- dim(inputTestData)
    inputTestLengthC    <- inputTestDataDim[1]
    inputTestLengthR    <- inputTestDataDim[2]
    testDistanceList    <- vector()
    classes             <- vector()
    testRandomization   <- sample(inputTestLengthC)
    testHit             <- 0
    testMiss            <- 0
    testGrand           <- 0
    print("=========================================")
    print("-- Classifying Data ---------------------")
    print("=========================================")
    for (i in testRandomization) {
        testGrand           <- testGrand + 1
        if (verboseClassify == TRUE) {
            print("-- New Sample ---------------------------")
            print("   Forward-Propagating...")
            print(paste("       Test Sample", i))
        }
        # Forward Propagation
        inputTestSample     <- inputTestData[i,]
        testHiddenLayer     <- sigmoid(inputTestSample %*% calibratedWeights)
        testHiddenLayer     <- c(insert(testHiddenLayer, element = 1.0)) # Insert Bias Term to Hidden Layer
        testOutput          <- sigmoid(testHiddenLayer %*% calibratedHidden)
        classes             <- c(classes, testOutput)
        # Metric Computation & Communication      
        if (verboseClassify == TRUE) print(paste("       Computed Class:        ", decimals(testOutput, 8)))
        if (length(testTargets > 0)) {
            testRounded         <- round(testOutput)
            testError           <- 0.5 * (testTargets[i] - testOutput) ^ 2
            testDistance        <- abs(testTargets[i] - testOutput) # An Easily Interpretable Error Measure
            testDistanceList    <- c(testDistanceList, testDistance) # For Monitoring Error Reduction  
            if (verboseClassify == TRUE) {
                print(paste("       Known Class:           ", testTargets[i]))
                print(paste("       Raw Distance:          ", decimals(testDistance, 8)))
                print(paste("       Computed SSE:          ", decimals(testError, 8)))
            }
            if (abs(testRounded - testTargets[i]) == 0) {
                testHit     <- testHit + 1
                if (verboseClassify == TRUE) print("       Rounded Hit Status:     Hit! :)")
            } else {
                testMiss    <- testMiss + 1
                if (verboseClassify == TRUE) print("       Rounded Hit Status:     Miss :(")
            }
            if (verboseClassify == TRUE) {
                print(paste("           Test Hits / Total: ", testHit, "/", testGrand))
                print(paste("           Test Hit Percent:  ", decimals((testHit/testGrand) * 100, 2)))
                print("-- Sample Done --------------------------")
            }
        }
    }
    print("=========================================")
    print("-- Classifications Complete -------------")
    print("=========================================")
    if (verboseClassify == FALSE) {
        if (length(testTargets > 0)) {
            print("-- Report -------------------------------")
            print("   Rounded Hits:")
            print(paste("       Test Hits / Total:     ", testHit, "/", testGrand))
            print(paste("       Test Hit Percent:      ", decimals((testHit/testGrand) * 100, 2)))
        }
    }
    # Return Results
    list (
        pValues     = classes,
        classes     = round(classes),
        hits        = testHit,
        total       = testGrand,
        percent     = decimals((testHit/testGrand) * 100, 2)
    )
}

# 
# Import Data
# 

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
testExpressionData.mod<-expressionData[which(!(names(expressionData) %in% testedCellLines.mod))]

#We need to reorder the columns according to trainKey
trainExpressionData.mod.reorder<-trainExpressionData.mod[, testedCellLines]
testExpressionData.mod.reorder<-testExpressionData.mod[, predictCellLines]

# Mash Data Frames Together & Remove Duplicate Row Names From trainKey 
trainingData        <- data.frame(cbind(trainKey[-1], t(as.data.frame(trainExpressionData.mod.reorder))))
# And Testing Data - Naming Conventions (Not Sure if This Will Actually be Needed - Does MLInterfaces handle train/validate/test?)
testingData         <- data.frame(t(as.data.frame(testExpressionData.mod.reorder)))

# Remove Row & Column Names & Convert to Matrix
names(trainExpressionData.mod.reorder)  <- NULL
names(testExpressionData.mod.reorder)   <- NULL
trainingData        <- t(data.matrix(trainExpressionData.mod.reorder))
testingData         <- t(data.matrix(testExpressionData.mod.reorder))

# 
# Setup Parameters
# 

# Standard Sample Selection
SampleSize          <- 40
Epochs              <- 50
EtaP                <- 0.08
EtaH                <- 0.2
HiddenNodes         <- 200
DWeightLimit        <- 0.5
HWeightLimit        <- 0.5

# High Variance Sample Selection
HighCV              <- TRUE
MinCV               <- 16
HiddenNodeRatio     <- 0.20
PlotTrainROC        <- FALSE

# 
# Additional Setup
# 

drugs               <- gsub("-", ".", trainKeyTransposed$Drug)
nDrugs              <- length(drugs)
tPredictions        <- list()
tSubset             <- vector()

# 
# Final Data Formatting
# 

# Create Training Subset & Add Bias Terms
trainingDim         <- dim(trainingData)
trainBiasTerms      <- generateVector(trainingDim[1], 1)
if (HighCV == TRUE) {
    tSubset         <- getByCV(trainingData, MinCV)
    hvTrainData     <- tSubset$matrix
    trainingSubset  <- cbind(trainBiasTerms, hvTrainData)
    trDim           <- dim(trainingSubset)
    HiddenNodes     <- round(trDim[2] * HiddenNodeRatio)
    SampleSize      <- trDim[2]
} else {
    trainingSubset  <- cbind(trainBiasTerms, trainingData[1:trainingDim[1],1:SampleSize])
}
colnames(trainingSubset) <- NULL

# Create Testing Subset & Add Bias Terms
testingDim          <- dim(testingData)
testBiasTerms       <- generateVector(testingDim[1], 1)
if (HighCV == TRUE) {
    hvTestData      <- testingData[ ,tSubset$columns]
    testingSubset   <- cbind(testBiasTerms, hvTestData)

} else {
    testingSubset   <- cbind(testBiasTerms, testingData[1:testingDim[1], 1:SampleSize])
}
colnames(testingSubset) <- NULL

# 
# Classifications & Predictions
# 

P <- proc.time()
for (i in 1:nDrugs) {
    print(paste("Drug ", i, ": ", drugs[i]), sep="")
    knownClasses    <- vector()
    for (known in trainKey[drugs[i]][,1]) knownClasses<- c(knownClasses, as.numeric(known))
    # Train Data
    computedValues  <- ANN.train (
        epochs              = Epochs,
        inputTrainingData   = trainingSubset,
        trainingTargets     = trainKey[drugs[i]][,1],
        etaP                = EtaP, 
        etaH                = EtaH, 
        hiddenNodes         = HiddenNodes, 
        dataWeightsLimit    = DWeightLimit, 
        hiddenWeightsLimit  = HWeightLimit, 
        plotData            = list(
            SSE                 = FALSE,
            distance            = FALSE,
            weightMeans         = list(
                plot                = FALSE,
                lables              = ""
            )
        ),
        verboseTrain        = FALSE
    )
    # Plot Training ROC
    if (PlotTrainROC == TRUE) {
        AUC             <- prediction(computedValues$classifications, knownClasses)
        AUCPerf         <- performance(AUC, "tpr", "fpr")
        plot(AUCPerf)
    }
    # Classify Test Data 
    tPredict        <- ANN.classify (
        inputTestData       = testingSubset, 
        calibratedWeights   = computedValues$trainedWeights, 
        calibratedHidden    = computedValues$trainedHidden,
        verboseClassify     = FALSE
    )
    # Test Predictions
    tPredictions    <- cbind(tPredictions, p = list(tPredict$classes))
}
RT           <- proc.time() - P
print(paste("Neural Net Runtime:", RT[3], "seconds"))

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
filename     <- paste("submissions/nn/genes_", SampleSize, "_epochs_", Epochs , "_hiddenNodes_", HiddenNodes, "_etaP_", EtaP, "_etaH_", EtaH, "_dataWeightsLimit_", DWeightLimit, "_hiddenWeightsLimit_", HWeightLimit, "_highCV_", HighCV, ".csv", sep="")
write.csv(submission, file=filename, row.names = FALSE)
print(paste("Output Saved As:", filename))

print("Status: Done")