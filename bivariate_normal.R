library(MASS)
library(MCMCpack)
library(graph)
library(ggplot2)
bluecrab <- read.table('/Users/jaideep/Documents/School/598/HW5/bluecrab.dat', header = F)
orangecrab <- read.table('/Users/jaideep/Documents/School/598/HW5/orangecrab.dat', header = F)
# Data for two different types of crabs
ybar.blue <- as.numeric(apply(bluecrab,2,mean))
ybar.orange <- as.numeric(apply(orangecrab,2,mean))

# Blue Crabs.
# Setup the Gibbs sampler
mu0 <- ybar.blue; S0 <- L0 <- sigma.bluecrab <- cov(bluecrab); n <- dim(bluecrab)[1]; nu0 <- 4;
THETA.BLUE <- SIGMA.BLUE <- NULL

set.seed(42)
for (s in 1:10000) {
  # update theta
  Ln <- solve(solve(L0) + n*solve(sigma.bluecrab))
  mun <- Ln%*%(solve(L0)%*%mu0 + n*solve(sigma.bluecrab)%*%ybar.blue)
  theta.bluecrab <- mvrnorm (1,mun,Ln)
  
  # update Sigma from above theta
  Sn <- S0 + ((t(bluecrab) - c(theta.bluecrab))%*%t(t(bluecrab) - c(theta.bluecrab)))
  sigma.bluecrab <- solve((rWishart(1,nu0+n,solve(Sn)))[,,1])
  
  THETA.BLUE <- rbind(THETA.BLUE, theta.bluecrab)
  SIGMA.BLUE <- rbind(SIGMA.BLUE, c(sigma.bluecrab))
  
}

plot(THETA.BLUE, main = 'Blue Crabs data', xlab = 'body depth', ylab = 'rear width')
abline(0,1)

# Covariance matracies created
rho.blue <- NULL
for (i in 1:dim(SIGMA.BLUE)[1]){
  rho.blue[i] <- cov2cor(matrix(SIGMA.BLUE[i,], 2, 2))[1,2]
}
plot(rho.blue, type = 'l')

# ORANGE CRABS

mu0 <- ybar.orange; S0 <- L0 <- sigma.orangecrab <- cov(orangecrab); n <- dim(orangecrab)[1]; nu0 <- 4;

THETA.ORANGE <- SIGMA.ORANGE <- NULL

set.seed(42)
for (s in 1:10000) {
  # update theta
  Ln <- solve(solve(L0) + n*solve(sigma.orangecrab))
  mun <- Ln%*%(solve(L0)%*%mu0 + n*solve(sigma.orangecrab)%*%ybar.orange)
  theta.orangecrab <- mvrnorm(1,mun,Ln)
  
  # update Sigma
  Sn <- S0 + ((t(orangecrab) - c(theta.orangecrab))%*%t(t(orangecrab) - c(theta.orangecrab)))
  sigma.orangecrab <- solve((rWishart(1,nu0+n,solve(Sn)))[,,1])
  
  THETA.ORANGE <- rbind(THETA.ORANGE, theta.orangecrab)
  SIGMA.ORANGE <- rbind(SIGMA.ORANGE, c(sigma.orangecrab))
  
}

plot(THETA.ORANGE, main = 'Orange Crabs data', xlab = 'body depth', ylab = 'rear width')
abline(0,1)


rho.orange <- NULL
for (i in 1:dim(SIGMA.ORANGE)[1]){
  rho.orange[i] <- cov2cor(matrix(SIGMA.ORANGE[i,], 2, 2))[1,2]
}

plot(rho.orange, type = 'l')
