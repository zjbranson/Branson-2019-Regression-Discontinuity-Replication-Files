library(MASS)
library(rstan)


rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#function that fits two Gaussian process regressions,
#one on the left side of the discontinuity, and one to the right.
#Requires the stan file twoGPsOnlyCovarianceParamsWithSlope.stan
#ARGUMENTS:
# data: a data.frame that has the running variable denoted as "x"
#       and an outcome denoted as "y"
# boundary: the value of x that denotes where the discontinuity is; default is 0.
# iters: number of iterations/draws for each chain in the MCMC for the Bayesian model
#        (which is run in rstan).
#        default is iters = 1000 for demonstration, but for actual application
#        this should be much higher, say 50,000 or 100,000.
# chains: number of chains to run in the MCMC for the Bayesian model. default is 2.
# Note: The total number of posterior draws this provides is (iters*chains)/2;
#       it's divided by 2 because half of the draws from each chain are thrown
#       out after burn-in.
fitTwoGPs = function(data, boundary = 0, iters = 1000, chains = 2){
  leftIndex = which(data$x < boundary)
  rightIndex = which(data$x >= boundary)
  x1 = data$x[leftIndex]; y1 = data$y[leftIndex]
  x2 = data$x[rightIndex]; y2 = data$y[rightIndex]

  gpFit = stan(file = "twoGPsOnlyCovarianceParamsWithSlope.stan",
                      data = list(x1 = x1, y1 = y1,
                        x2 = x2, y2 = y2,
                        N1 = length(x1), N2 = length(x2)),
                      iter= iters,
                      chains = chains)
    return(gpFit)
}
#After the Bayesian model using fitTwoGPs() is run,
#this function performs kriging (or extrapolation) to obtain
#posterior draws for the treatment effect at the boundary in an RDD.
# n: for each draw of the posterior, the number of draws from the predictive
#    posterior distribution of the Gaussian process that are obtained. default is 1,
#    which is fine for inference.
# data: a data.frame that has the running variable denoted as "x"
#       and an outcome denoted as "y"
# stanFit: the object returned by fitTwoGPs()
# boundary: the value of x that denotes where the discontinuity is; default is 0.
# length: the number of points specified for kriging on either side of the boundary.
#         "length" does not affect the precision of inference for the average treatment effect,
#         but making length larger can be useful for visualizations of the mean response functions.
#         "length" is only used in the getPredX1() and getPredX2() functions;
#         more details about these functions are below.
# Note: this function returns a matrix with (n*iters*chains/2)-many rows and
#       (Nt + Nc + 2*length)-many columns,
#where Nt and Nc are the number of units in treatment and control, respectively.
# Importantly, the (Nc + length)-th column denotes the posterior draws for the response
# at the boundary on the left-hand side (assuming the left-hand side is the control),
# and the (Nc + length + 1)-th column denotes the posterior draws for the response
# at the boundary on the right-hand side.
performKrigingSameParams = function(n = 1, data, stanFit, boundary = 0, length = 10){
  x1 = data$x[which(data$x < boundary)]
  y1 = data$y[which(data$x < boundary)]
  x2 = data$x[which(data$x >= boundary)]
  y2 = data$y[which(data$x >= boundary)]
  
  predX1 = getPredX1(data = data, boundary = boundary, length = length)
  predX2 = getPredX2(data = data, boundary = boundary, length = length)

  #extract the variables from the stanFit object
  stanFit = extract(stanFit)
  
  #first, perform kriging on the left side
  krigingLeft = drawPostPredY(n = n,
                mu = stanFit$mu1,
                beta = stanFit$beta1,
                  sigmasq = stanFit$sigmasq,
                              phi = stanFit$phi,
                              etasq = stanFit$etasq,
                              X = x1,
                              Y = y1,
                              predX = predX1)
  #now perform kriging on the right side
  krigingRight = drawPostPredY(n = n,
                  mu = stanFit$mu2,
                  beta = stanFit$beta2,
                  sigmasq = stanFit$sigmasq,
                              phi = stanFit$phi,
                              etasq = stanFit$etasq,
                              X = x2,
                              Y = y2,
                              predX = predX2)
  
  #krigingLeft and krigingRight are matrices
  #We'll cbind these matrices, so that
  #The first 100 columns are for the left,
  #the next 100 columns for the right.
  krigingMatrix = cbind(krigingLeft, krigingRight)
  return(krigingMatrix)
}
#draw from the posterior predictive distribution of the Gaussian process.
#Details about this process can be found in the Branson et al. RDD paper.
drawPostPredY = function(n = 1, mu, sigmasq, phi, etasq, X, Y, predX,beta){
  
  #First we need distances among the X values
  distances = matrix(nrow=length(X), ncol=length(X))
  for(i in 1:length(X)){
    for(j in 1:length(X)){
      distances[i,j] = X[i] - X[j]
    }
  }
  #and now we need the distances among the *new* X values
  predDistances = matrix(nrow=length(predX), ncol=length(predX))
  for(i in 1:length(predX)){
    for(j in 1:length(predX)){
      predDistances[i,j] = predX[i] - predX[j]
    }
  }
  #and also we need the distances between each of the X values and the *new* X values
  predXDistances = t(sapply(X, function(x) x - predX))
  
  #sigmasq and etasq are paramters whose posterior we already drew from
  #using the Stan model.
  #Thus, for each posterior draw of the parameters, we'll draw from the posterior
  #predictive distribution for each new X (i.e., each predX)
  #Thus, we'll have length(sigmasq)-many draws from the posterior predictive distribution.
  posteriorPredictiveDraws = list(length = length(sigmasq))
  for(m in 1:length(sigmasq) ){
    print(m)
    
    #the covariance matrix for the Xs is
    estCovMat = sigmasq[m]*exp(-phi[m]*distances^2) + diag(etasq[m], length(X))
    #and the covariance matrix for the new Xs is
    predCovMat = sigmasq[m]*exp(-phi[m]*predDistances^2)
    #and the covariance matrix between the Xs and new Xs is
    predCrossCovMat = sigmasq[m]*exp(-phi[m]*predXDistances^2)
    
    #using conditional MVN theory, we can find the distribution of
    #p(predX | X, Y, theta)
    #the mean of this distribution is
    predYMean = (mu[m] + predX*beta[m]) + t(predCrossCovMat)%*%solve(estCovMat)%*%(Y - mu[m] - X*beta[m] )
    #and the covariance matrix of this distribution is
    predYCovMat = predCovMat - t(predCrossCovMat)%*%solve(estCovMat)%*%predCrossCovMat
    #Therefore, using the above mean and covariance, we can draw from the posterior
    #predictive distribution for each predX.
    posteriorPredictiveDraws[[m]] = mvrnorm(n = n, mu = predYMean, Sigma = predYCovMat)
  }
  #right now, posteriorPredictiveDraws is a list of matrices,
  #and we'd like to collapse this into one big matrix
  posteriorPredictiveDrawsMatrix = do.call(rbind, posteriorPredictiveDraws)
  return(posteriorPredictiveDrawsMatrix)
} 

#function that obtains data points to be used for kriging
#Essentially, this function chooses "length" many points between
#the right-most point on the left-hand-side of the boundary,
#as well as "length" many points between
#the left-most point on the right-hand-side of the boundary.
#This then returns Nt + length + Nc + length points,
#where Nt and Nc are the number of units in treatment and control, respectively.
getPredX1 = function(data, boundary, length = 10){
  x1 = data$x[which(data$x < boundary)]
  predX1 = c(x1, seq(max(x1), boundary, length = length))
  return(predX1)
}
getPredX2 = function(data, boundary, length = 10){
  x2 = data$x[which(data$x >= boundary)]
  predX2 = c(seq(boundary, min(x2), length = length), x2)
  return(predX2)
}

pickMeans.mp = readRDS("pickMeans.mp.rds")

#in this example, the boundary is 30.5
boundary = 30.5
#the number of treatment and control units is
Nc = length(which(pickMeans.mp$x <= boundary))
Nt = length(which(pickMeans.mp$x > boundary))
length = 10

#obtain the Bayesian model fit from rstan
fit.test = fitTwoGPs(data = pickMeans.mp, iters = 1000, boundary = boundary)
#perform kriging (i.e., extrapolation of the Gaussian process regressions to the boundary)
kriging.test = performKrigingSameParams(n = 1, data = pickMeans.mp, stanFit = fit.test, length = length, boundary = boundary)
#finally, obtain posterior draws of the treatment effect.
#This just corresponds to the difference in posterior draws for
# the response at the boundary on the right-hand side and
# the response at the boundary on the left-hand side
treatmentEffect.test = kriging.test[,Nc + length + 1] - kriging.test[,Nc + length]
#the corresponding point estimate and confidence interval should be close to
#what is reported for Minutes Played in Table 1 in the paper:
#(note that in the paper we set iters = 50000, which will take a bit more time to run.)
mean(treatmentEffect.test)
quantile(treatmentEffect.test, probs = c(0.025, 0.975))