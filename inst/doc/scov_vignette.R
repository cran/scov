## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = FALSE)
library(mvtnorm)
library(scov)

## ----intercept----------------------------------------------------------------
intercept = matrix(1,ncol=4,nrow=4)
corrplot::corrplot(intercept,method = "square")
mtext("intercept", side = 2, line = 0, cex = 1)

## ----X1-----------------------------------------------------------------------
X1 = rbind(c(1,1,1,0),c(1,1,1,0),c(1,1,1,0),c(0,0,0,1))
corrplot::corrplot(X1,method = "square")
mtext("X1", side = 2, line = 0, cex = 1)

## ----X2-----------------------------------------------------------------------
X2 = rbind(c(1,0,0,0),c(0,1,1,1),c(0,1,1,1),c(0,1,1,1))
corrplot::corrplot(X2,method = "square")
mtext("X2", side = 2, line = 0, cex = 1)

## ----interaction--------------------------------------------------------------
corrplot::corrplot(X2*X1,method = "square")
mtext("interaction between X1 and X2", side = 2, line = 0, cex = 1)

## ----covariate list-----------------------------------------------------------
covar_mats = list(intercept=intercept,X1=X1,X2=X2)

## ----spatial------------------------------------------------------------------
adj_matrix = rbind(c(0,1,0,0),c(1,0,0,0),c(0,0,0,1),c(0,0,1,0))
corrplot::corrplot(adj_matrix,method = "square")
mtext("adjacency matrix", side = 2, line = 0, cex = 1)

## ----load data----------------------------------------------------------------
mean = rep(0,4)
sigma = 0.05*intercept+0.2*X1+0.2*X2+0.1*X2*X1+0.4*(diag(4) + adj_matrix)
diag(sigma) = 1
num_observations = 100
dataset = mvtnorm::rmvnorm(num_observations,mean=mean,sigma=sigma)

## ----show sigma---------------------------------------------------------------
corrplot::corrplot(sigma,method = "square")
mtext("correlation matrix", side = 2, line = 0, cex = 1)

## ----compute sce--------------------------------------------------------------
sce = scov::scov(covar_mats, dataset, adj_matrix,
                 interaction_effects=list(c("X1","X2")),
                 ncores=1,parallelize=FALSE,verbose=FALSE)
corrplot::corrplot(sce$corrmat_estim,method = "square")
mtext("SCE", side = 2, line = 0, cex = 1)

## ----initialize non-normal data-----------------------------------------------
dataset = mvtnorm::rmvt(num_observations,sigma=sigma,df=2) + matrix(runif(4*num_observations,max=10001,min=10000),nrow=num_observations,ncol=4)

## ----compute ive--------------------------------------------------------------
ive = scov::scov(covar_mats, dataset, adj_matrix,
                 interaction_effects=list(c("X1","X2")),
                 semiparametric=TRUE,
                 ncores=1,parallelize=FALSE)
corrplot::corrplot(ive$corrmat_estim,method = "square")
mtext("IVE", side = 2, line = 0, cex = 1)

## ----load data 2--------------------------------------------------------------
mean = rep(0,4)
sigma = 0.05*intercept+0.2*X1+0.2*X2+0.1*X2*X1+0.4*(diag(4) + adj_matrix)
diag(sigma) = 1
dataset = mvtnorm::rmvnorm(num_observations,mean=mean,sigma=sigma)

## ----compute wsce-------------------------------------------------------------
covar_mats = list(intercept=intercept,X1=X1)
wsce = scov::scov(list(X1=X1), dataset, adj_matrix,
                  misspecification = TRUE,
                  parallelize = FALSE,
                  ncores=1,verbose=FALSE)
corrplot::corrplot(wsce$corrmat_estim,method = "square")
mtext("WSCE", side = 2, line = 0, cex = 1)

## ----show lambda--------------------------------------------------------------
paste0("lambda: ",wsce$lambda)

