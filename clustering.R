#Clustering script
#Written by Jason Li

#hierarchical clustering
d=dist(as.matrix(t(expressionData[,-1])))
hc<-hclust(d)
plot(hc, main="Cluster Dendrogram of Whole Expression Set")
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
kmfit <- kmeans(t(expressionData[,-1]), 4)
kmfitsub<-kmeans(t(expressionData.sub), 4)
kmfitpca<-kmeans(expr.pca$x, 4)

# Cluster Plot against 1st 2 principal components

# vary parameters for most readable graph
library(cluster) 
#pammy<-pam(t(expressionData[,-1]), 4, metric="euclidean")
clusplot(t(expressionData[,-1]), kmfit$cluster, color=TRUE, shade=TRUE, labels=2, lines=0)

# Centroid Plot against 1st 2 discriminant functions
library(fpc)
#plotcluster(t(expressionData[,-1]), kmfit$cluster) #VERY SLOW, do not run!

expr.pca<-prcomp(t(expressionData[,-1])) #pca of whole ES
exprsub.pca<-prcomp(t(expressionData.sub)) #pca of whole ES + subtype
#plot kmeans 
plot(expr.pca$x[,1], expr.pca$x[,2], col=(as.integer(kmfit$cluster)), pch=20, 
     xlab="PC1", ylab="PC2", main="K-means clustering")
text(expr.pca$x[,1], expr.pca$x[,2], names(expressionData[,-1]), col=as.integer(kmfit$cluster))


#plot(exprsub.pca$x[,1], exprsub.pca$x[,2], col=(as.integer(kmfitsub$cluster)), pch=20, 
#     xlab="PCA1", ylab="PCA2", main="K-means clustering")

#plot(expr.pca$x[,3], expr.pca$x[,4], col=(as.integer(expressionData.sub[18633,])), pch=20, 
#     xlab="PC1", ylab="PC2", main="Subtype")
#Plot PC, color by Subtype
plot(exprsub.pca$x[,1], exprsub.pca$x[,2], col=(as.integer(expressionData.sub[18633,])), pch=20, 
     xlab="PC1", ylab="PC2", main="Subtype")
text(expr.pca$x[,1], expr.pca$x[,2], names(expressionData[,-1]), col=as.integer(kmfit$cluster))

#Plot PC, color by kmeans
plot(exprsub.pca$x[,1], exprsub.pca$x[,2], col=(as.integer(expressionData.sub[18633,])), pch=20, 
     xlab="PC1", ylab="PC2", main="Subtype")
text(exprsub.pca$x[,1], exprsub.pca$x[,2], names(expressionData.sub), col=as.integer(expressionData.sub[18633,]))

plot(exprsub.pca$x[,1], exprsub.pca$x[,2], col=(as.integer(kmfitsub$cluster)), pch=20, 
     xlab="PC1", ylab="PC2", main="K-means clustering with subtype data")
text(exprsub.pca$x[,1], exprsub.pca$x[,2], names(expressionData.sub), col=as.integer(kmfitsub$cluster))
