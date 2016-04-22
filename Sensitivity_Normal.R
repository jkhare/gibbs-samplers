# Sensitivity analysis.
# Prior parameters.
mu0.sens <- 75; s20.sens <- 100; n <- 16;
# vector of prior hyper-parameters
prec.prior <- matrix(c(1,2,4,8,16,32,64,128,256,1,2,4,8,16,32,64,128,256), ncol = 2)

# data properties + initialize vector of probabilities
ybar.a <- 75.2; ybar.b <- 77.5; s.a <- 7.3; s.b <- 8.1;
pr <- NULL

# looping through different values of hyper-parameters of the prior k0 and nu0 and building a Gibbs sampler for each model
set.seed(42)
for (i in 1:dim(prec.prior)[1]){
  k0 <- prec.prior[i,1]; nu0 <- prec.prior[i,2];
  kn <- k0 + n; nun <- nu0 + n;
  mun.a <- (k0*mu0.sens + n*ybar.a)/kn
  s2n.a <- (nu0*s20.sens + (n-1)*(s.a^2) + k0*n*(ybar.a-mu0.sens)^2/(kn))/nun
  s2.post.a <- 1/rgamma(10000, nun/2, s2n.a*nun/2)
  theta.post.a <- rnorm(1000, mun.a, sqrt(s2.post.a/kn))
  mun.b <- (k0*mu0.sens + n*ybar.b)/kn
  s2n.b <- (nu0*s20.sens + (n-1)*(s.b^2) + k0*n*((ybar.b-mu0.sens)^2/(kn))/nun)
  s2.post.b <- 1/rgamma(10000, nun/2, s2n.b*nun/2)
  theta.post.b <- rnorm(1000, mun.b, sqrt(s2.post.b/kn))
  # Find the probability that the posterior mean of a is less than that of b
  pr <- c(pr, mean(theta.post.a<theta.post.b))
}

# vector of probabilities that the mean of sample a is less than mean of sample b, plotted against the range of values of hyper-parameters
plot(prec.prior[,1],pr,type = 'b', main = 'probabilities v/s (k0 = nu0)', xlab = 'values of (k0=nu0)', ylab = 'probabilities')
