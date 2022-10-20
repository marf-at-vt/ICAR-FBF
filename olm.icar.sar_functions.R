#### OLM, ICAR, and SAR fractional integrated likelihoods
#### Erica M. Porter

# fractional integrated likelihood for OLM
# this integrated likelihood can be found analytically so does not require
# adaptive quadrature or other techniques
vague.logmargin.independent.q <- function(X, Y, q=0.1) {

  p <- ncol(X)
  n <- length(Y)

  logmarginal.q <- (0.5*(n*(q-1)))*log(2*pi) + (0.5*p)*log(q) + lgamma(0.5*(n-p)) -
    lgamma(0.5*(n*q-p)) + (0.5*(p-n))*log(0.5*(t(Y)%*%(diag(n) - X%*%solve(t(X)%*%X)%*%t(X))%*%Y)) +
    (0.5*(n*q-p))*log(0.5*(q*(t(Y)%*%(diag(n) - X%*%solve(t(X)%*%X)%*%t(X))%*%Y)))

  return(logmarginal.q)
}

############################## SAR model with independence Jeffreys prior ##############################

# integrand over tau_c for the fractional integrated likelihood for SAR model
# requires adjacency matrix W of 0's and 1's indicating if each pair of regions are neighbors
sar.integrand.spatial.q <- function(phi, q=1, Y, X, W, logmax=0) {

  result <- rep(NA,length(phi))
  integrand <- rep(NA,length(phi))
  n <- length(Y)
  num.reg <- length(Y)
  p <- ncol(X)
  norm.constant <- ((0.5*(p-n*q))*log(2*pi)) - ((0.5*p)*log(q)) + (lgamma((0.5*(n*q-p))))

  # eigenvalues of W
  lambda.eig <- eigen(W, symmetric=T)$values

  # Set up info with Sig.phi.inverse matrix
  for(i in 1:length(phi)){
    Sig.phi.inverse <- (diag(num.reg) - phi[i]*W)%*%(diag(num.reg) - phi[i]*W)
    s2.phi <- t(Y)%*%(Sig.phi.inverse - Sig.phi.inverse%*%X%*%solve(t(X)%*%Sig.phi.inverse%*%X)%*%t(X)%*%Sig.phi.inverse)%*%Y

    # set up prior for phi
    sum1 <- sum((lambda.eig/(1-phi[i]*lambda.eig))^2)
    sum2 <- sum((lambda.eig/(1-phi[i]*lambda.eig)))
    phi.prior.log <- 0.5*log((sum1 - (1/n)*(sum2^2)))

    # Fractional integrand (with all known normalizing constants)

    integrand[i] <- norm.constant + (q/2)*determinant(Sig.phi.inverse,logarithm=TRUE)$modulus -
      0.5*log(det(t(X)%*%Sig.phi.inverse%*%X)) +
      ((p-n*q)/2)*log((q/2)*s2.phi) + phi.prior.log - logmax

    result[i] <- exp(integrand[i])
  }
  return(result)
}

# function to help avoid numerical underflow during integration
sar.logmax.spatial.q <- function(phi, q=1, Y, X, W, logmax=0) {

  integrand <- rep(NA,length(phi))
  n <- length(Y)
  num.reg <- length(Y)
  p <- ncol(X)
  norm.constant <- ((0.5*(p-n*q))*log(2*pi)) - ((0.5*p)*log(q)) + (lgamma((0.5*(n*q-p))))

  # eigenvalues of W
  lambda.eig <- eigen(W, symmetric=T)$values

  # Set up info with Sig.phi.inverse matrix
  for(i in 1:length(phi)){
    Sig.phi.inverse <- (diag(num.reg) - phi[i]*W)%*%(diag(num.reg) - phi[i]*W)
    s2.phi <- t(Y)%*%(Sig.phi.inverse - Sig.phi.inverse%*%X%*%solve(t(X)%*%Sig.phi.inverse%*%X)%*%t(X)%*%Sig.phi.inverse)%*%Y

    # set up prior for phi
    sum1 <- sum((lambda.eig/(1-phi[i]*lambda.eig))^2)
    sum2 <- sum((lambda.eig/(1-phi[i]*lambda.eig)))
    phi.prior.log <- 0.5*log((sum1 - (1/n)*(sum2^2)))

    # Fractional integrand (with all known normalizing constants)

    integrand[i] <- norm.constant + (q/2)*determinant(Sig.phi.inverse,logarithm=TRUE)$modulus -
      0.5*log(det(t(X)%*%Sig.phi.inverse%*%X)) +
      ((p-n*q)/2)*log((q/2)*s2.phi) + phi.prior.log - logmax
  }
  return(integrand)
}

############################## ICAR model with reference prior ##############################

# integrand over tau_c for the fractional integrated likelihood for ICAR model
# requires neighborhood matrix H
ref.integrand.spatial.q <- function(tauc,q=0.1,Y,X,H,H.spectral=NULL,Sig_phi=NULL,logmax) {

  result <- rep(NA,length(tauc))
  n <- length(Y)
  num.reg <- length(Y)
  p <- ncol(X)
  norm.constant <- (0.5*(p-n*q))*log(2*pi) - (0.5*p)*log(q)

  # Get info from adjacency matrix H
  if(is.null(H.spectral)==TRUE){
    row.names(H) <- NULL
    H.spectral <- eigen(H,symmetric=TRUE)
  }

  Q <- H.spectral$vectors
  eigH <- H.spectral$values

  if(is.null(Sig_phi)==TRUE){
    phimat <- diag(1/sqrt(eigH[1:(num.reg-1)]))
    Sig_phi <- matrix(0,num.reg,num.reg) #initialize
    for(i in 1:(num.reg-1)){
      total <- (1/(eigH[i]))*Q[,i]%*%t(Q[,i])
      Sig_phi <- Sig_phi + total
    }
  }

  # Projection of X and spectral decomp
  s1 <- diag(num.reg) - X%*%solve(t(X)%*%X)%*%t(X)
  s1 <- (s1+t(s1))/2
  s1.eigen <- eigen(s1, symmetric = TRUE)
  Q.star <- s1.eigen$vectors[,1:(num.reg-ncol(X))]
  #Q.star <- eigen(s1)$vectors[,1:(num.reg-ncol(X))]
  M <- t(Q.star)%*%Sig_phi%*%Q.star
  M.spectral <- eigen(M, symmetric = TRUE)
  U <- M.spectral$vectors
  xi <- M.spectral$values
  L <- Q.star%*%U

  # Get chi-square quantity for approximating S2 (line 16 Bayes paper)
  LtY.squared <- (t(L)%*%Y)^2

  #Y.sp <- t(Q) %*% Y
  X.sp <- t(Q) %*% X

  for(i in 1:length(tauc)){

    # Set up info for omega matrix
    aux1 <- c((1+(1/(tauc[i]*(eigH[1:(n-1)])))),1)
    #omega.det <- prod(aux1)
    omega.log.det <- sum(log(aux1))
    Ainv <- matrix(1/aux1, nrow=n, ncol=ncol(X))
    s2.tauc <- sum(LtY.squared/(1+((tauc[i]^(-1))*xi)))

    aux0 <- t(X.sp*Ainv)
    #aux <- det(aux0 %*% X.sp)
    aux <- determinant(aux0 %*% X.sp,logarithm=T)$modulus[1]
    aux2 <- xi/(tauc[i]+xi)

    # Fractional integrand (with all normalizing constants)
    integrand <- norm.constant - (0.5*q)*omega.log.det -
      (0.5*aux) - log(tauc[i]) +
      lgamma((0.5*(n*q-p))) + (0.5*(p-n*q))*log(0.5*(q*s2.tauc)) +
      (0.5*log((sum((aux2)^2) - (1/(n-p))*((sum((aux2)))^2)))) - logmax

    result[i] <- exp(integrand)
  }
  return(result)
}

# function to help avoid numerical underflow during integration
ref.logmax.spatial.q <- function(tauc,q=0.1,Y,X,H,H.spectral=NULL,Sig_phi=NULL) {

  integrand <- rep(NA,length(tauc))
  n <- length(Y)
  num.reg <- length(Y)
  p <- ncol(X)
  norm.constant <- (0.5*(p-n*q))*log(2*pi) - (0.5*p)*log(q)

  # Get info from adjacency matrix H
  if(is.null(H.spectral)==TRUE){
    row.names(H) <- NULL
    H.spectral <- eigen(H,symmetric=TRUE)
  }

  Q <- H.spectral$vectors
  eigH <- H.spectral$values

  if(is.null(Sig_phi)==TRUE){
    phimat <- diag(1/sqrt(eigH[1:(num.reg-1)]))
    Sig_phi <- matrix(0,num.reg,num.reg) #initialize
    for(i in 1:(num.reg-1)){
      total <- (1/(eigH[i]))*Q[,i]%*%t(Q[,i])
      Sig_phi <- Sig_phi + total
    }
  }

  # Projection of X and spectral decomp
  s1 <- diag(num.reg) - X%*%solve(t(X)%*%X)%*%t(X)
  s1 <- (s1+t(s1))/2
  s1.eigen <- eigen(s1, symmetric = TRUE)
  Q.star <- s1.eigen$vectors[,1:(num.reg-ncol(X))]
  #Q.star <- eigen(s1)$vectors[,1:(num.reg-ncol(X))]
  M <- t(Q.star)%*%Sig_phi%*%Q.star
  M.spectral <- eigen(M, symmetric = TRUE)
  U <- M.spectral$vectors
  xi <- M.spectral$values
  L <- Q.star%*%U

  # Get chi-square quantity for approximating S2 (line 16 Bayes paper)
  LtY.squared <- (t(L)%*%Y)^2

  #Y.sp <- t(Q) %*% Y
  X.sp <- t(Q) %*% X

  # Set up info for omega matrix
  for(i in 1:length(tauc)){
    aux1 <- c((1+(1/(tauc[i]*(eigH[1:(n-1)])))),1)
    #omega.det <- prod(aux1)
    omega.log.det <- sum(log(aux1))
    Ainv <- matrix(1/aux1, nrow=n, ncol=ncol(X))
    s2.tauc <- sum(LtY.squared/(1+((tauc[i]^(-1))*xi)))

    aux0 <- t(X.sp*Ainv)
    #aux <- det(aux0 %*% X.sp)
    #aux <- det(aux0 %*% X.sp,log=T)
    aux <- determinant(aux0 %*% X.sp,logarithm=T)$modulus[1]
    aux2 <- xi/(tauc[i]+xi)

    # Fractional integrand (with all normalizing constants)
    integrand[i] <- norm.constant - (0.5*q)*omega.log.det -
      (0.5*aux) - log(tauc[i]) +
      lgamma((0.5*(n*q-p))) + (0.5*(p-n*q))*log(0.5*(q*s2.tauc)) +
      (0.5*log((sum((aux2)^2) - (1/(n-p))*((sum((aux2)))^2))))
  }
  return(integrand)
}

############################## FBF model selection ##############################
############################## with all 3 model types ###########################

probs.sar.icar <- function(Y,X,W,H,H.spectral=NULL,Sig_phi=NULL,q=0.05) {

  X <- as.matrix(X)
  n <- length(Y)
  p <- ncol(X)-1
  num.mods <- (2^p)*3
  lambda.1.inv <- 1/eigen(W,symmetric = T)$values[1]
  lambda.n.inv <- 1/eigen(W,symmetric = T)$values[n]
  phi.grid <- seq(lambda.n.inv+0.01,lambda.1.inv-0.01,by=.001)
  tauc.grid <- c(seq(0.001,10,by=.01),seq(10.1,100,by=10))

  ind.logmargin <- rep(NA,(num.mods/3))
  icar.logmargin <- rep(NA,(num.mods/3))
  sar.logmargin <- rep(NA,(num.mods/3))

  icar.logmargin.1 <- rep(NA,(num.mods/3))
  icar.logmargin.q <- rep(NA,(num.mods/3))
  sar.logmargin.1 <- rep(NA,(num.mods/3))
  sar.logmargin.q <- rep(NA,(num.mods/3))
  ind.BF <- rep(NA,(num.mods/3))
  icar.BF <- rep(NA,(num.mods/3))
  sar.BF <- rep(NA,(num.mods/3))

  ind.mod.probs <- rep(NA,(num.mods/3))
  icar.mod.probs <- rep(NA,(num.mods/3))
  sar.mod.probs <- rep(NA,(num.mods/3))

  if(p>1){
    ## dX is a matrix of 0's and 1's indicating all possible covariate combinations
    x <- rbind(rep(0,p),combos(p)$binary)
    dX <- cbind(rep(1,2^p),x)

    ## Create a vector of formula labels for each candidate model
    predictor <- lapply(1:2^p, matrix, data=NA)
    predictor.labels <- lapply(1:2^p,matrix, data=NA)
    predictor.formulas <- lapply(1:2^p,matrix, data=NA)
    formula.vec <- c(rep(NA, 2^p))
    labels <- c(rep(NA,p))
    for(i in 1:(p)){labels[i] <- paste("X",i,sep="")}
    labels <- c("Intercept",labels)

    for(i in 1:2^p){
      predictor[[i]] <- diag(dX[i,])
      colnames(predictor[[i]]) <- labels
      predictor.labels[[i]]<- c(rownames(as.matrix((which(colSums(predictor[[i]]) != 0)))))
      predictor.formulas[[i]] <- formula(paste("Y ~ ", paste(predictor.labels[[i]], collapse=" + ")))
      formula.vec[i] <- Reduce(paste,deparse(predictor.formulas[[i]]))
    }

    mod.prior <- rep(NA,2^p)
    for(i in 1:2^p){
      k <- (length(which(colSums(predictor[[i]]) != 0)) - 1)
      mod.prior[i] <- ((1/(p+1))*(nchoosek(p,k))^(-1))/2
    }

    mod.prior <- c(mod.prior,mod.prior/2,mod.prior/2)

    ## combo.list is a list of matrices with the actual covariate values for each combination in dX
    ## Consider the case of p>1, p=1, and p=0 separately because the combos function accommodates p>1
    combo.list <- lapply(1:2^p, matrix, data= NA)
    for(i in 1:2^p){
      combo.list[[i]] <- X%*%diag(dX[i,])
    }
    for(i in 1:((2^p)-1)) {combo.list[[i]] <- as.matrix(combo.list[[i]][,-(which(colSums(combo.list[[i]]) == 0))])}
  }

  else if(p==1) {
    formula.vec <- c(rep(NA, 2^p))
    formula.vec[1] <- Reduce(paste,deparse(formula(paste("Y ~ ", paste("Intercept", collapse=" + ")))))
    formula.vec[2] <- Reduce(paste,deparse(formula(paste("Y ~ ", paste("Intercept + X1", collapse=" + ")))))

    combo.list <- lapply(1:2^p, matrix, data= NA)
    combo.list[[2]] <- X
    combo.list[[1]] <- as.matrix(X[,1])

    mod.prior <- rep(NA,2^p)
    for(i in 1:2^p){
      k <- 1
      mod.prior[i] <- ((1/(p+1))*(nchoosek(p,k))^(-1))/2
    }

    mod.prior <- c(mod.prior,mod.prior/2,mod.prior/2)
  } else if(p==0) {
    formula.vec <- c(rep(NA, 2^p))
    formula.vec[1] <- Reduce(paste,deparse(formula(paste("Y ~ ", paste("Intercept", collapse=" + ")))))

    combo.list <- lapply(1:2^p, matrix, data= NA)
    combo.list[[1]] <- X
    mod.prior <- c(0.5,0.25,0.25)
  }

  icar.joint <- array(dim=c(2^p,2,length(tauc.grid)))
  sar.joint <- array(dim=c(2^p,2,length(phi.grid)))
  icar.logmax.value <- array(dim=c(2^p,2))
  sar.logmax.value <- array(dim=c(2^p,2))

  for(i in 1:(num.mods/3)){
    ind.logmargin[i] <- vague.logmargin.independent.q(q=q, X=combo.list[[i]], Y=Y)

    sar.joint[i,1,] <- sar.logmax.spatial.q(phi.grid,Y=Y,X=combo.list[[i]],W=W,q=1)

    sar.joint[i,2,] <- sar.logmax.spatial.q(phi.grid,Y=Y,X=combo.list[[i]],W=W,q=q)

    sar.logmax.value[i,1] <- max(na.omit(sar.joint[i,1,]))
    sar.logmax.value[i,2] <- max(na.omit(sar.joint[i,2,]))

    sar.logmargin.1[i] <- log(integrate(sar.integrand.spatial.q,lower=lambda.n.inv,upper=lambda.1.inv,
                                        Y=Y,X=combo.list[[i]],W=W,
                                        q=1,logmax=sar.logmax.value[i,1])$value) + sar.logmax.value[i,1]

    sar.logmargin.q[i] <- log(integrate(sar.integrand.spatial.q,lower=lambda.n.inv,upper=lambda.1.inv,
                                        Y=Y,X=combo.list[[i]],W=W,
                                        q=q,logmax=sar.logmax.value[i,2])$value) + sar.logmax.value[i,2]

    sar.logmargin[i] <- sar.logmargin.1[i] - sar.logmargin.q[i]

    print(paste('SAR Integrated Likelihood ', i,' of ', num.mods/3, ' complete',' at ',date(),sep=''))
  }

  for(i in 1:(num.mods/3)){
    icar.joint[i,1,] <- ref.logmax.spatial.q(tauc.grid,Y=Y,X=combo.list[[i]],H=H,
                                             H.spectral=H.spectral,Sig_phi=Sig_phi,q=1)

    icar.joint[i,2,] <- ref.logmax.spatial.q(tauc.grid,Y=Y,X=combo.list[[i]],H=H,
                                             H.spectral=H.spectral,Sig_phi=Sig_phi,q=q)

    icar.logmax.value[i,1] <- max(icar.joint[i,1,])
    icar.logmax.value[i,2] <- max(icar.joint[i,2,])

    icar.logmargin.1[i] <- log(integrate(ref.integrand.spatial.q,lower=0,upper=Inf,
                                         Y=Y,X=combo.list[[i]],H=H,H.spectral=H.spectral,
                                         Sig_phi=Sig_phi,q=1,logmax=icar.logmax.value[i,1])$value) + icar.logmax.value[i,1]

    icar.logmargin.q[i] <- log(integrate(ref.integrand.spatial.q,lower=0,upper=Inf,
                                         Y=Y,X=combo.list[[i]],H=H,H.spectral=H.spectral,
                                         Sig_phi=Sig_phi,q=q,logmax=icar.logmax.value[i,2])$value) + icar.logmax.value[i,2]

    icar.logmargin[i] <- icar.logmargin.1[i] - icar.logmargin.q[i]

    print(paste('ICAR Integrated Likelihood ', i,' of ', num.mods/3, ' complete',' at ',date(),sep=''))
  }

  logmargin.all <- c(ind.logmargin,sar.logmargin,icar.logmargin)
  base.model <- max(logmargin.all)

  for(i in 1:(num.mods/3)){
    ind.BF[i] <- exp(ind.logmargin[i]-base.model)
    sar.BF[i] <- exp(sar.logmargin[i]-base.model)
    icar.BF[i] <- exp(icar.logmargin[i]-base.model)
  }

  BF.vec <- c(ind.BF,sar.BF,icar.BF)
  BF.sum.adj <- BF.vec %*% mod.prior

  mod.probs.all <- rep(NA,num.mods)
  for (i in 1:num.mods){
    mod.probs.all[i] <- (BF.sum.adj)^(-1)*BF.vec[i]*mod.prior[i]
  }

  probs.mat <- data.frame(mod.probs.all,c(rep("Independent",num.mods/3),rep("SAR",num.mods/3),rep("ICAR",num.mods/3)),c(formula.vec,formula.vec,formula.vec))
  names(probs.mat) <- c("model prob","model type","model form")

  return(list(probs.mat=probs.mat,
              mod.prior=mod.prior,
              logmargin.all=logmargin.all,
              base.model=base.model,
              BF.vec=BF.vec))
}

# FBF selection using uniform model priors
probs.sar.icar_uniform <- function(Y,X,W,H,H.spectral=NULL,Sig_phi=NULL,q=0.05) {

  X <- as.matrix(X)
  n <- length(Y)
  p <- ncol(X)-1
  num.mods <- (2^p)*3
  lambda.1.inv <- 1/eigen(W,symmetric = T)$values[1]
  lambda.n.inv <- 1/eigen(W,symmetric = T)$values[n]
  phi.grid <- seq(lambda.n.inv+0.01,lambda.1.inv-0.01,by=.001)
  tauc.grid <- c(seq(0.001,10,by=.01),seq(10.1,100,by=10))

  ind.logmargin <- rep(NA,(num.mods/3))
  icar.logmargin <- rep(NA,(num.mods/3))
  sar.logmargin <- rep(NA,(num.mods/3))

  icar.logmargin.1 <- rep(NA,(num.mods/3))
  icar.logmargin.q <- rep(NA,(num.mods/3))
  sar.logmargin.1 <- rep(NA,(num.mods/3))
  sar.logmargin.q <- rep(NA,(num.mods/3))
  ind.BF <- rep(NA,(num.mods/3))
  icar.BF <- rep(NA,(num.mods/3))
  sar.BF <- rep(NA,(num.mods/3))

  ind.mod.probs <- rep(NA,(num.mods/3))
  icar.mod.probs <- rep(NA,(num.mods/3))
  sar.mod.probs <- rep(NA,(num.mods/3))

  if(p>1){
    ## dX is a matrix of 0's and 1's indicating all possible covariate combinations
    x <- rbind(rep(0,p),combos(p)$binary)
    dX <- cbind(rep(1,2^p),x)

    ## Create a vector of formula labels for each candidate model
    predictor <- lapply(1:2^p, matrix, data=NA)
    predictor.labels <- lapply(1:2^p,matrix, data=NA)
    predictor.formulas <- lapply(1:2^p,matrix, data=NA)
    formula.vec <- c(rep(NA, 2^p))
    labels <- c(rep(NA,p))
    for(i in 1:(p)){labels[i] <- paste("X",i,sep="")}
    labels <- c("Intercept",labels)

    for(i in 1:2^p){
      predictor[[i]] <- diag(dX[i,])
      colnames(predictor[[i]]) <- labels
      predictor.labels[[i]]<- c(rownames(as.matrix((which(colSums(predictor[[i]]) != 0)))))
      predictor.formulas[[i]] <- formula(paste("Y ~ ", paste(predictor.labels[[i]], collapse=" + ")))
      formula.vec[i] <- Reduce(paste,deparse(predictor.formulas[[i]]))
    }

    mod.prior <- rep(1/num.mods,num.mods)

    ## combo.list is a list of matrices with the actual covariate values for each combination in dX
    ## Consider the case of p>1, p=1, and p=0 separately because the combos function accommodates p>1
    combo.list <- lapply(1:2^p, matrix, data= NA)
    for(i in 1:2^p){
      combo.list[[i]] <- X%*%diag(dX[i,])
    }
    for(i in 1:((2^p)-1)) {combo.list[[i]] <- as.matrix(combo.list[[i]][,-(which(colSums(combo.list[[i]]) == 0))])}
  }

  else if(p==1) {
    formula.vec <- c(rep(NA, 2^p))
    formula.vec[1] <- Reduce(paste,deparse(formula(paste("Y ~ ", paste("Intercept", collapse=" + ")))))
    formula.vec[2] <- Reduce(paste,deparse(formula(paste("Y ~ ", paste("Intercept + X1", collapse=" + ")))))

    combo.list <- lapply(1:2^p, matrix, data= NA)
    combo.list[[2]] <- X
    combo.list[[1]] <- as.matrix(X[,1])

    mod.prior <- rep(1/num.mods,num.mods)

  } else if(p==0) {
    formula.vec <- c(rep(NA, 2^p))
    formula.vec[1] <- Reduce(paste,deparse(formula(paste("Y ~ ", paste("Intercept", collapse=" + ")))))

    combo.list <- lapply(1:2^p, matrix, data= NA)
    combo.list[[1]] <- X
    mod.prior <- c(1/3,1/3,1/3)
  }

  icar.joint <- array(dim=c(2^p,2,length(tauc.grid)))
  sar.joint <- array(dim=c(2^p,2,length(phi.grid)))
  icar.logmax.value <- array(dim=c(2^p,2))
  sar.logmax.value <- array(dim=c(2^p,2))

  for(i in 1:(num.mods/3)){
    ind.logmargin[i] <- vague.logmargin.independent.q(q=q, X=combo.list[[i]], Y=Y, a1=0, b1=0)

    sar.joint[i,1,] <- sar.logmax.spatial.q(phi.grid,Y=Y,X=combo.list[[i]],W=W,q=1)

    sar.joint[i,2,] <- sar.logmax.spatial.q(phi.grid,Y=Y,X=combo.list[[i]],W=W,q=q)

    sar.logmax.value[i,1] <- max(na.omit(sar.joint[i,1,]))
    sar.logmax.value[i,2] <- max(na.omit(sar.joint[i,2,]))

    sar.logmargin.1[i] <- log(integrate(sar.integrand.spatial.q,lower=lambda.n.inv,upper=lambda.1.inv,
                                        Y=Y,X=combo.list[[i]],W=W,
                                        q=1,logmax=sar.logmax.value[i,1])$value) + sar.logmax.value[i,1]

    sar.logmargin.q[i] <- log(integrate(sar.integrand.spatial.q,lower=lambda.n.inv,upper=lambda.1.inv,
                                        Y=Y,X=combo.list[[i]],W=W,
                                        q=q,logmax=sar.logmax.value[i,2])$value) + sar.logmax.value[i,2]

    sar.logmargin[i] <- sar.logmargin.1[i] - sar.logmargin.q[i]

    print(paste('SAR Integrated Likelihood ', i,' of ', num.mods/3, ' complete',' at ',date(),sep=''))
  }

  for(i in 1:(num.mods/3)){
    icar.joint[i,1,] <- ref.logmax.spatial.q(tauc.grid,Y=Y,X=combo.list[[i]],H=H,
                                             H.spectral=H.spectral,Sig_phi=Sig_phi,q=1)

    icar.joint[i,2,] <- ref.logmax.spatial.q(tauc.grid,Y=Y,X=combo.list[[i]],H=H,
                                             H.spectral=H.spectral,Sig_phi=Sig_phi,q=q)

    icar.logmax.value[i,1] <- max(icar.joint[i,1,])
    icar.logmax.value[i,2] <- max(icar.joint[i,2,])

    icar.logmargin.1[i] <- log(integrate(ref.integrand.spatial.q,lower=0,upper=Inf,
                                         Y=Y,X=combo.list[[i]],H=H,H.spectral=H.spectral,
                                         Sig_phi=Sig_phi,q=1,logmax=icar.logmax.value[i,1])$value) + icar.logmax.value[i,1]

    icar.logmargin.q[i] <- log(integrate(ref.integrand.spatial.q,lower=0,upper=Inf,
                                         Y=Y,X=combo.list[[i]],H=H,H.spectral=H.spectral,
                                         Sig_phi=Sig_phi,q=q,logmax=icar.logmax.value[i,2])$value) + icar.logmax.value[i,2]

    icar.logmargin[i] <- icar.logmargin.1[i] - icar.logmargin.q[i]

    print(paste('Ref Integrated Likelihood ', i,' of ', num.mods/3, ' complete',' at ',date(),sep=''))
  }


  logmargin.all <- c(ind.logmargin,sar.logmargin,icar.logmargin)
  base.model <- max(logmargin.all)

  for(i in 1:(num.mods/3)){
    ind.BF[i] <- exp(ind.logmargin[i]-base.model)
    sar.BF[i] <- exp(sar.logmargin[i]-base.model)
    icar.BF[i] <- exp(icar.logmargin[i]-base.model)
  }

  BF.vec <- c(ind.BF,sar.BF,icar.BF)
  BF.sum.adj <- BF.vec %*% mod.prior

  mod.probs.all <- rep(NA,num.mods)
  for (i in 1:num.mods){
    mod.probs.all[i] <- (BF.sum.adj)^(-1)*BF.vec[i]*mod.prior[i]
  }

  probs.mat <- data.frame(mod.probs.all,c(rep("Independent",num.mods/3),rep("SAR",num.mods/3),rep("ICAR",num.mods/3)),c(formula.vec,formula.vec,formula.vec))
  names(probs.mat) <- c("model prob","model type","model form")

  return(list(probs.mat=probs.mat,
              mod.prior=mod.prior,
              logmargin.all=logmargin.all,
              base.model=base.model,
              BF.vec=BF.vec))
}

