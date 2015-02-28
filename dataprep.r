# Load Libraries 
# source("http://bioconductor.org/biocLite.R")
# biocLite("MLInterfaces")
# install.packages("MLInterfaces")
library(MLInterfaces)

# Load Data
subtypeData         <- read.csv("data/subtype.csv")
idMapData           <- read.csv("data/scoring_and_test_set_id_mappings.csv")
# Can Select on DFs Below Using DF$ColumnName List Syntax 
trainKey            <- read.csv("data/training_set_answers_transposed.csv")
expressionData      <- read.csv("data/expression.csv")

# Set Variables 
cellLineNames       <- names(trainKey)

# Loops to Merge Expression Data & Training Key Together
for (i in 2:ncol(expressionData)) {
    for (j in 2:ncol(trainKey)) {
        if (cellLineNames[j] == names(expressionData)[i]) {
            
        }
        #     print(head(trainKey[,j]))
    }
#     print(head(expressionData[,i]))
}
# Outer join: merge(x = df1, y = df2, by = "CustomerId", all = TRUE)
# 
# Left outer: merge(x = df1, y = df2, by = "CustomerId", all.x=TRUE)
# 
# Right outer: merge(x = df1, y = df2, by = "CustomerId", all.y=TRUE)
# 
# Cross join: merge(x = df1, y = df2, by = NULL)