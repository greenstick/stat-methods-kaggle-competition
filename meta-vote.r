setwd("~/Documents/Education/Graduate/OHSU/Courses/Winter 2015/Statistical Methods/assignments/project1/src") 

# Crudely gathers files from every nn & svm run and does a majority vote for response of each cell line
nnFiles <- Sys.glob("submissions/nn/*.csv")
svmFiles<- Sys.glob("submissions/svm/*.csv")
files 	<- c(nnFiles, svmFiles)
aFile 	<- read.csv(files[1])$value
n 		<- length(files)
l 		<- length(aFile)
data 	<- c(rep(0, l))
for (file in files) {
	d <- read.csv(file)$value
	data <- data + d
}
≤f
ids          <- seq(1, l, 1)
values 		 <- round(data / n)
subdf        <- data.frame(id = as.matrix(ids), value = as.matrix(values))
submission   <- data.frame(lapply(subdf, as.character), stringsAsFactors=FALSE)
filename     <- paste("submissions/meta/", n, "files.csv", sep="")
write.csv(submission, file=filename, row.names = FALSE)
print(paste("Output Saved As:", filename))