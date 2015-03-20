setwd("~/Documents/Education/Graduate/OHSU/Courses/Winter 2015/Statistical Methods/assignments/project1/src") 

# Crude Bootstrapping Function
bootstrap 	 <- function (data, sampleSize = 30) {
	sample  	<- sample(1:dim(data)[2], sampleSize)
	out 		<- vector()
	for (value in sample) {
		d 		<- data[,value]
		out 	<- cbind(out, d)
	}
	out 	 <- rowMeans(out)
	out
}

# Crude Bagging Function
bagData 	<- function (data, samplesize, iterations) {
	out 	<- vector()
	for (i in 1:iterations) {
		boot 	<- bootstrap(data = data, sampleSize = samplesize)
		out  	<- cbind(out, boot)
	}
	out <- rowMeans(out)
	out
}

# Set Some Params
highScoring  <- FALSE
ss 			 <- 3
iter 		 <- 2000

# Crudely gathers files from every nn & svm run and does a majority vote for response of each cell line
if (highScoring == TRUE) {
	files 	<- Sys.glob("submissions/high-score/*.csv")
} else {
	nnFiles <- Sys.glob("submissions/nn/*.csv")
	svmFiles<- Sys.glob("submissions/svm/*.csv")
	files 	<- c(nnFiles, svmFiles)
}
aFile 	<- read.csv(files[1])$value
n 		<- length(files)
l 		<- length(aFile)
data 	<- vector()
for (file in files) {
	d <- read.csv(file)$value
	if (length(d) == l) {
		data <- cbind(data, d)
	}
}

# Pretty Crude, Bootstrap Samples all SVM & Neural Net Classification Files, then Bags the Results of the voted classification
baggery 	 <- bagData(data = data, samplesize = 120, iterations = 1000)


ids          <- seq(1, l, 1)
values 	     <- baggery
rValues      <- round(baggery)
subdf        <- data.frame(id = as.matrix(ids), value = as.matrix(values))
rSubdf       <- data.frame(id = as.matrix(ids), value = as.matrix(rValues))
submission   <- data.frame(lapply(subdf, as.character), stringsAsFactors=FALSE)
rSubmission  <- data.frame(lapply(rSubdf, as.character), stringsAsFactors=FALSE)
filename     <- paste("submissions/meta/", n, "-iterations", iter, "-samplesize", ss, "files.csv", sep="")
rFilename 	 <- paste("submissions/meta/", n, "-iterations", iter, "-samplesize", ss, "files-rounded.csv", sep="")
write.csv(submission, file=filename, row.names = FALSE)
write.csv(rSubmission, file=rFilename, row.names = FALSE)
print(paste("Output Saved As:", filename))
print(paste("Output Saved As:", rFilename))