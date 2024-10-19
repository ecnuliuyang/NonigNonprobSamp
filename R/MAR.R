##' Doubly robust inference with non-probability survey samples
##'
##' @description Function to apply the pseudo-likelihood method in the assumption of a ignorable participation mechanism.
##'
##' @param SA A matrix, the non-probability sample. The last column indicates the outcome variable and the other columns represent the covariates.
##' @param SB A matrix, the reference probability sample. The last column indicates the survey weights and the other columns represent the covariates.
##'
##' @return A list with two elements:
##' \itemize{
##'   \item mu, the IPW, regression, and double robust estimates of the population mean.
##'   \item theta, the maximum pseudo-likelihood estimates of parameters in the participation probability model.
##' }
##'
##'
##' @references
##' Chen, Y., Li, P. and Wu, C. (2020).
##' Doubly robust inference with nonprobability survey samples.
##' \emph{Journal of the American Statistical Association} \strong{115}, 2011–2021.
##'
##' @importFrom stats plogis glm binomial
##'
##' @export
##'
MAR = function(SA,SB){


  xA = SA[,-ncol(SA)]
  yA = SA[,ncol(SA)]
  xB = SB[,-ncol(SB)]
  d = SB[,ncol(SB)]

  mcoef = glm(y~. , data = data.frame(SA),family = binomial(link = "logit"))$coefficients
  hatyA = as.numeric(plogis(cbind(1,xA)%*%mcoef))
  hatyB = as.numeric(plogis(cbind(1,xB)%*%mcoef))

  eta = rep(0,ncol(xB)+1)
  err = 1
  while (err>0.00001) {
    deltaB = as.numeric(cbind(1,xB)%*%eta)
    pi = exp(deltaB)/(1+exp(deltaB))
    g = colSums(cbind(1,xA)) - colSums(d*pi*cbind(1,xB))
    h = t(cbind(1,xB))%*%(d*pi*(1-pi)*cbind(1,xB))
    upd = solve(h)%*%g
    eta= eta+upd
    err = max(abs(upd))
  }
  theta_hat= as.numeric(eta)

  piA_hat = as.numeric(exp(cbind(1,xA)%*%theta_hat)/(1+exp(cbind(1,xA)%*%theta_hat)))

  #IPW
  hatNA = sum(1/piA_hat)
  mu1 = sum(yA/piA_hat)/hatNA
  #Regression
  mu2 = sum(d*hatyB)/sum(d)
  #DR
  mu3 = sum((yA -hatyA)/piA_hat)/hatNA + mu2

  list(mu = c(mu1,mu2,mu3), theta = theta_hat)
}
