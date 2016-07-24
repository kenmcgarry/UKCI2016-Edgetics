# UKCI2016-reviewers.r
# Sharon McDonald
# random graphs and ROC curves for UKCI-2016
# The two reviewers felt that further confidence could be placed on our networks if we compared 
# them with random graph models. They also suggested using ROC curves are used to demonstrate accuracy.

library(igraph)
library(xtable)
library(sand)
library(ROCR)
library(eigenmodel)

setwd("C:/R-files/UKCI2016")
# ensure that PPI is in environment by loading UKCI2016-reviewers.RData

gr1 <- sample_gnm(225, 2594, directed = FALSE, loops = FALSE) # 225 proteins with 2594 connections between them
gr1 <- simplify(gr1)

nedges<-2594
nverts<- 225
avepath <- average.path.length(gr1)
connected<- is.connected(gr1)
avedeg <- mean(degree(gr1))
diam<- diameter(gr1,weights=NA)
modular<- modularity(gr1, membership(cluster_walktrap(gr1)))
transit=transitivity(gr1)

randomg <- data.frame(modular,avepath,nedges,nverts,transit,avedeg,diam,connected)

xtable(randomg) # the single line of results was pasted into the original latex table No 2


# ----- ROC curves from Kolaczyk book (page 107-108) ----
# ensure that PPI is in environment by loading workspace "UKCI2016-reviewers.RData"

# CHUNK 1
# Now need to bind number of snps, number of interactions and if a disease gene to dataframe
colnames(ppi) <- c("protein", "partner","confidence")
y <- as.data.frame(nacts)
x <- as.data.frame(nsnps)
z <- merge(x,y,by="protein")
mu <- merge(ppi,z,by="protein")

PPI <-graph.data.frame(ppi,directed = FALSE)
A <- get.adjacency(PPI,sparse=FALSE)
dim(A)

v.attrs <- get.data.frame(PPI, what="vertices")
# CHUNK 30
perm.index <- sample(1:17020) # 185 * 184/2= 17112 permutations
nfolds <- 5
nmiss <- 17020/nfolds

Avec <- A[lower.tri(A)]
Avec.pred1 <- numeric(length(Avec))

# CHUNK 31
for(i in seq(1,nfolds)){
  # Index of missing values.
  miss.index <- seq(((i-1) * nmiss + 1),(i*nmiss), 1)
  A.miss.index <- perm.index[miss.index]
  
  # Fill a new Atemp appropriately with NA's.
  Avec.temp <- Avec
  Avec.temp[A.miss.index] <- rep("NA", length(A.miss.index))
  Avec.temp <- as.numeric(Avec.temp)
  Atemp <- matrix(0, 185, 185) # varying these two numbers leads to different accuracies
  Atemp[lower.tri(Atemp)] <- Avec.temp
  Atemp <- Atemp + t(Atemp)
  
  # Now fit model and predict.
  Y <- Atemp
  model1.fit <- eigenmodel_mcmc(Y, R=2,S=11000, burn=10000)
  model1.pred <- model1.fit$Y_postmean
  model1.pred.vec <- model1.pred[lower.tri(model1.pred)]
  Avec.pred1[A.miss.index] <- model1.pred.vec[A.miss.index]
}

# CHUNK 32
pred1 <- prediction(Avec.pred1, Avec)
perf1 <- performance(pred1, "tpr", "fpr")

plot(perfall, col="red", lwd=5)
plot(perf1, add = TRUE, col="blue",lwd=5)
abline(a=0, b= 1,lty=2 )

# CHUNK 33
perf1.auc <- performance(pred1, "auc")
slot(perf1.auc, "y.values")





