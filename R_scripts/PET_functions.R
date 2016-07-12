############################################
# Pleiotropy Estimation and Testing (PET)  #
# Functions by Zhang et al Genet Epidemiol #
############################################

estimateTraitResidual <- function(data, trait, 
                                  covars, link, distribution) {
  # calculate the covariate adjusted residuals
  # from a linear model for a trait
  
  t.covars <- paste(covars, collapse=" + ", sep=" + ")
  t.form <- paste(trait, t.covars, sep=" ~ ")
  t.mod <- glm(t.form, family=distribution(link), data=data)
  
  t.residuals <- t.mod$residuals
  
  return(t.residuals)  
}

constructDataMatrix <- function(trait1, trait2, geno) {
  # create the matrix for the LMM with X, T and IT components
  
  n = length(geno)
  design.vec <- c(rep(1, n), rep(2, n))
  Y <- c(trait1, trait2)
  X <- c(geno, geno)
  
  I.vec <- design.vec * X
  mu <- mean(trait1 + trait2)
  
  input.matrix <- data.frame(cbind(1, Y, X, design.vec, I.vec))
  
  colnames(input.matrix) <- c("mu", "Y", "X", "Tvec", "Ivec")
  
  return(input.matrix)
}


calculateResidualCovariance <- function(data, method="ML"){
  # calculate residual covariance from LMM residuals
  require(nlme)
  lmm.fit <- lme(Y ~ X + Tvec + Ivec, random=~1|mu, data=data,
                 na.action=na.omit,
                 method=method)
  
  # the returned residuals is a 2-segment vector
  covar.mat <- cov(matrix(resid(lmm.fit), ncol=2))
  return(covar.mat[1, 2])
}


calcPleiotropyCorrelation <- function(trait1, trait2, residual.covar){
  # calculate the pleiotropy correlation coefficient (PCC) from
  # the difference between the observed trait covariance and
  # residual model covariance, standardized to trait standard deviations
  
  trait.covar <- cov(trait1, trait2)
  trait1.sd <- sd(trait1)
  trait2.sd <- sd(trait2)
  
  delta <- trait.covar - residual.covar
  rho <- abs((delta)/(trait1.sd * trait2.sd))

  return(rho)  
}

PleiotropyEstimationTest <- function(data.set, trait1, trait2,
                                     genotypes){
  
  # calculate the pleiotropy correlation coefficient
  t1 <- data.set[[trait1]]
  t2 <- data.set[[trait2]]
  geno <- data.set[[genotypes]]
  
  input.matrix <- constructDataMatrix(t1, t2, geno)
  resid.covar <- calculateResidualCovariance(input.matrix, "ML")
  
  pcc <- calcPleiotropyCorrelation(t1, t2, resid.covar)
  
  return(pcc)
}

BootstrapPET <- function(data, indices, trait1, trait2, genotypes){
  data.mat <- data[indices, ]
  pcc = PleiotropyEstimationTest(data.mat, trait1, trait2, genotypes)
  pcc
}

PETB <- function(data, trait1, trait2, genotypes, resamples,
                 plot=FALSE){
  # calculate the assymetric 2-sided P-value from bootstrap resampling
  require(boot)
  
  boot.out <- boot(data=data, statistic=BootstrapPET, R=resamples,
                   trait1=trait1, trait2=trait2, genotypes=genotypes)
  if(plot){
    plot(boot.out)
  }
  
  pval <- mean(abs(boot.out$t0) > abs(boot.out$t))
  
  return(c(boot.out$t0, pval))
}