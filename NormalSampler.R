# Hierarchical modeling of school data
school <- NULL
school[1] <- read.table('/Users/jaideep/Documents/School/598/HW4/school1.dat', header = F)
school[2] <- read.table('/Users/jaideep/Documents/School/598/HW4/school2.dat', header = F)
school[3] <- read.table('/Users/jaideep/Documents/School/598/HW4/school3.dat', header = F)
school[4] <- read.table('/Users/jaideep/Documents/School/598/HW4/school4.dat', header = F)
school[5] <- read.table('/Users/jaideep/Documents/School/598/HW4/school5.dat', header = F)
school[6] <- read.table('/Users/jaideep/Documents/School/598/HW4/school6.dat', header = F)
school[7] <- read.table('/Users/jaideep/Documents/School/598/HW4/school7.dat', header = F)
school[8] <- read.table('/Users/jaideep/Documents/School/598/HW4/school8.dat', header = F)

# Data parameters
ybar.school <- sapply(school, mean)
m <- length(school)
var.school <- sapply(school, var)
n.school <- sapply(school, length)

# prior parameters
mu0.school <- 7; g20.school <- 5; tau20.school <- 10; eta0.school <- 2; s20.school <- 15; nu0.school <- 2;

# Gibbs Sampler
# Initialization
GIBBSLENGTH <- 10000
set.seed(42)
THETA.FINAL.SCHOOL <- matrix(nrow = S, ncol = m)
SMT.FINAL.SCHOOL <- matrix(nrow = S, ncol = 3)
sigma2.school <- mean(var.school)
mu.school <- mean(ybar.school)
tau2.school <- var(ybar.school)
theta.school <- mean(ybar.school)

# Sample values
for (s in 1:GIBBSLENGTH){
  # first sample theta
  for (j in 1:m) {
    vtheta.school <- 1/(n.school[j]/sigma2.school + 1/tau2.school)
    etheta.school <- vtheta.school*(ybar.school[j]*n.school[j]/sigma2.school + mu.school/tau2.school)
    theta.school[j] <- rnorm(1,etheta.school,sqrt(vtheta.school))
  }
  # update sigma2
  nun.school <- nu0.school + sum(n.school)
  ss.school <- nu0.school*s20.school
  for (j in 1:m){
    ss.school <- ss.school + sum((as.data.frame(school[j]) - theta.school[j])^2)
  }
  sigma2.school <- 1/rgamma(1,nun.school/2,ss.school/2)

  #update mu
  vmu.school <- 1/((m/tau2.school) + 1/g20.school)
  emu.school <- vmu.school*((m*mean(theta.school)/tau2.school) + (mu0.school/g20.school))
  mu.school <- rnorm (1,emu.school,sqrt(vmu.school))

  # update tau2
  etam.school <- eta0.school + m
  ss.school <- (eta0.school*tau20.school + sum ((theta.school - mu.school)^2))
  tau2.school <- 1/rgamma(1,etam.school/2,ss.school/2)

  # store generated values
  THETA.FINAL.SCHOOL[s,] <- theta.school
  SMT.FINAL.SCHOOL[s,] <- c(sigma2.school,mu.school,tau2.school)
}

# Compare sample means with posterior means
plot(colMeans(THETA.FINAL.SCHOOL), ybar.school)
# Compare sample mean with hyper-parameter
mean(ybar.school)
mean(SMT.FINAL.SCHOOL[,2])

# Posterior density of R
R<-NULL
for (i in 1:dim(SMT.FINAL.SCHOOL)[1]){
  R[i] <- SMT.FINAL.SCHOOL[i,3]/(SMT.FINAL.SCHOOL[i,1] + SMT.FINAL.SCHOOL[i,3])
}
plot(R)
