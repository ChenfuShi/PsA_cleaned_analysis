library("apeglm")
library("emdbook")
options(warn=-1)

# Load the data
cts <- read.csv(".local/hic_goingtoR/hic_CD4_tot_counts.csv", header=TRUE, row.names=1)
ase.cts <- read.csv(".local/hic_goingtoR/hic_CD4_ASE_counts.csv", header=TRUE, row.names=1)
cts <- as.matrix(cts)
ase.cts <- as.matrix(ase.cts)
mean_values = rowMeans(cts, na.rm = TRUE)

# fill NA
weights <- ifelse(is.na(ase.cts), 1e-6, 1)
ase.cts[is.na(ase.cts)] <- 1
cts[is.na(cts)] <- 2

theta.hat <- 100 # rough initial estimate of dispersion
X=matrix(rep(1,ncol(cts)))
niter <- 3
system.time({
  for (i in 1:niter) {
    param <- cbind(theta.hat, cts)
    fit.mle <- apeglm(Y=ase.cts, x=X, log.lik=NULL, param=param,
                      no.shrink=TRUE, log.link=FALSE, method="betabinCR", weights = weights)
    print(i)
    print("fit mle done")
    theta.hat <- bbEstDisp(success=ase.cts, size=cts,
                           x=X, beta=fit.mle$map,
                           minDisp=1, maxDisp=500, weights = weights)
    print(i)
    print("fit dispersion done")
  }
})

coef <- 1
mle <- cbind(fit.mle$map[,coef], fit.mle$sd[,coef])
param <- cbind(theta.hat, cts)
system.time({
  fit2 <- apeglm(Y=ase.cts, x=X, log.lik=NULL, param=param,
                 coef=coef, mle=mle, threshold=0.5,
                 log.link=FALSE, method="betabinCR", weights = weights)
})
s.val <- svalue(fit2$thresh)

# save to text file
write.table(s.val, file=".local/hic_goingtoR/apeglm_results_CD4.txt", sep=",")




####################################################################################################




# Load the data
cts <- read.csv(".local/hic_goingtoR/hic_CD8_tot_counts.csv", header=TRUE, row.names=1)
ase.cts <- read.csv(".local/hic_goingtoR/hic_CD8_ASE_counts.csv", header=TRUE, row.names=1)
cts <- as.matrix(cts)
ase.cts <- as.matrix(ase.cts)
mean_values = rowMeans(cts, na.rm = TRUE)

# fill NA
weights <- ifelse(is.na(ase.cts), 1e-6, 1)
ase.cts[is.na(ase.cts)] <- 1
cts[is.na(cts)] <- 2

theta.hat <- 100 # rough initial estimate of dispersion
X=matrix(rep(1,ncol(cts)))
niter <- 3
system.time({
  for (i in 1:niter) {
    param <- cbind(theta.hat, cts)
    fit.mle <- apeglm(Y=ase.cts, x=X, log.lik=NULL, param=param,
                      no.shrink=TRUE, log.link=FALSE, method="betabinCR", weights = weights)
    print(i)
    print("fit mle done")
    theta.hat <- bbEstDisp(success=ase.cts, size=cts,
                           x=X, beta=fit.mle$map,
                           minDisp=1, maxDisp=500, weights = weights)
    print(i)
    print("fit dispersion done")
  }
})

coef <- 1
mle <- cbind(fit.mle$map[,coef], fit.mle$sd[,coef])
param <- cbind(theta.hat, cts)
system.time({
  fit2 <- apeglm(Y=ase.cts, x=X, log.lik=NULL, param=param,
                 coef=coef, mle=mle, threshold=0.5,
                 log.link=FALSE, method="betabinCR", weights = weights)
})
s.val <- svalue(fit2$thresh)

# save to text file
write.table(s.val, file=".local/hic_goingtoR/apeglm_results_CD8.txt", sep=",")




####################################################################################################




# Load the data
cts <- read.csv(".local/hic_goingtoR/hic_ALL_tot_counts.csv", header=TRUE, row.names=1)
ase.cts <- read.csv(".local/hic_goingtoR/hic_ALL_ASE_counts.csv", header=TRUE, row.names=1)
cts <- as.matrix(cts)
ase.cts <- as.matrix(ase.cts)
mean_values = rowMeans(cts, na.rm = TRUE)

# fill NA
weights <- ifelse(is.na(ase.cts), 1e-6, 1)
ase.cts[is.na(ase.cts)] <- 1
cts[is.na(cts)] <- 2

theta.hat <- 100 # rough initial estimate of dispersion
X=matrix(rep(1,ncol(cts)))
niter <- 3
system.time({
  for (i in 1:niter) {
    param <- cbind(theta.hat, cts)
    fit.mle <- apeglm(Y=ase.cts, x=X, log.lik=NULL, param=param,
                      no.shrink=TRUE, log.link=FALSE, method="betabinCR", weights = weights)
    print(i)
    print("fit mle done")
    theta.hat <- bbEstDisp(success=ase.cts, size=cts,
                           x=X, beta=fit.mle$map,
                           minDisp=1, maxDisp=500, weights = weights)
    print(i)
    print("fit dispersion done")
  }
})

coef <- 1
mle <- cbind(fit.mle$map[,coef], fit.mle$sd[,coef])
param <- cbind(theta.hat, cts)
system.time({
  fit2 <- apeglm(Y=ase.cts, x=X, log.lik=NULL, param=param,
                 coef=coef, mle=mle, threshold=0.5,
                 log.link=FALSE, method="betabinCR", weights = weights)
})
s.val <- svalue(fit2$thresh)

# save to text file
write.table(s.val, file=".local/hic_goingtoR/apeglm_results_ALL.txt", sep=",")


