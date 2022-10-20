## FBF case study - Columbus, OH data set
## Response: crime rate

################################### Instructions #####################################
## Before running the code below, either copy and paste the functions from the file ##
## olm.icar.sar_functions.R into the R console or source the the R script from      ##
## a directory on your computer.                                                    ##
######################################################################################

# remove all environment variables
rm(list = ls())

# load R libraries for loading data and performing model selection
library(spdep)
library(sf)
library(hier.part)
library(pracma)

# read in the data as contained in the spdep package
columbus <- st_read(system.file("shapes/columbus.shp", package="spData")[1], quiet=TRUE)

# create neighborhood matrix
columbus.listw <- poly2nb(columbus)
summary(columbus.listw)
W <- nb2mat(columbus.listw, style="B")
Dmat <- diag(apply(W,1,sum))
num.reg <- length(columbus$CRIME)

H <- Dmat - W
H <- (H+t(H))/2
rownames(H) <- NULL
isSymmetric(H)

# spectral quantities for use in model selection
H.spectral <- eigen(H, symmetric=TRUE)
Q <- H.spectral$vectors
eigH <- H.spectral$values
phimat <- diag(1/sqrt(eigH[1:(num.reg-1)]))
Sig_phi <- matrix(0,num.reg, num.reg) #initialize
for(i in 1:(num.reg-1)){
  total <- (1/(eigH[i]))*Q[,i]%*%t(Q[,i])
  Sig_phi <- Sig_phi + total
}

# define response and design matrix
Y <- columbus$CRIME
X <- cbind(1, columbus$HOVAL, columbus$INC, columbus$OPEN, columbus$PLUMB, columbus$DISCBD)
q <- (ncol(X)+1)/num.reg

# source selection functions
source("olm.icar.sar_functions.R")

# run model selection for all candidate OLM, ICAR, and SAR models
columbus.select <- probs.sar.icar(Y=Y,X=X,W=W,H=H,
                                  H.spectral=H.spectral,
                                  Sig_phi=Sig_phi,
                                  q=q)

# print the model with highest posterior model probability
columbus.select$probs.mat[which.max(columbus.select$probs.mat[,1]),]

# print top 6 models
columbus.select$probs.mat[order(columbus.select$probs.mat$`model prob`, decreasing=TRUE),][1:6,]

# print vector of posterior inclusion probabilities covariates
post.include.cov <- matrix(NA,nrow = 1, ncol=ncol(X)-1)
labels <- c(rep(NA, ncol(X)-1))
for(i in 1:(ncol(X)-1)){labels[i] <- paste("X", i, sep="")}
colnames(post.include.cov) <- labels

for(j in 1:ncol(X)-1){
  post.include.cov[,j] <- sum(columbus.select$probs.mat[grep(paste("X",j,sep=""), columbus.select$probs.mat$'model form'), 1])
}

post.include.cov
