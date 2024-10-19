##' Doubly robust inference with nonprobability survey samples
##'
##' @description Function to apply the pseudo-likelihood method in the assumption of a nonignorable participation mechanism.
##'
##' @param SA A matrix, the non-probability  sample. The last column indicates the outcome variable and the other columns represent the covariates. Additionally, the second-to-last column serves as the instrumental variable.
##' @param SB A matrix, the reference probability sample. The last column indicates the survey weights and the other columns represent the covariates. Additionally, the second-to-last column serves as the instrumental variable.
##' @param par0 A vector, the initial value of the parameter vector in the participation probability model. The default value is set to zero.
##'
##' @return A list with three elements:
##' \itemize{
##'   \item theta, the maximum pseudo-likelihood estimates of the parameters in the participation probability model.
##'   \item xi, the maximum likelihood estimates of the parameters in the outcome regression model.
##'   \item mu, the regression, IPW, and AIPW estimates of the population mean.
##' }
##'
##'
##' @references
##' Liu, Y., Yuan, M., Li, P. and Wu, C. (2024).
##' Statistical inference with nonignorable non-probability survey samples.
##' \emph{arXiv}: 2410.02920.
##'
##' @importFrom stats plogis glm nlminb binomial
##'
##' @export
##'
MNAR <- function (SA, SB, par0 = rep(0,ncol(SA))) {

  ell <- function (theta) {

    alpha <- theta[1]
    beta <- theta[-c(1, length(theta))]
    gamma <- theta[length(theta)]

    qA <- c(uA%*%beta)
    qB <- c(uB%*%beta)
    muA <- c(xi1 + xA%*%xi2)
    cA <- log((1 + exp(muA+gamma))/(1 + exp(muA)))
    muB <- c(xi1 + xB%*%xi2)
    cB <- log((1 + exp(muB+gamma))/(1 + exp(muB)))

    sum(alpha + qA + cA) - sum(dB * (alpha + qB + cB)) +
      sum(dB * log(1 + exp(alpha + qB + cB)))
  }

  yA <- c(SA[,ncol(SA)])
  xA <- SA[, -ncol(SA)]
  uA <- xA[,-ncol(xA)]
  xB <- SB[,-ncol(SB)]
  dB <- SB[,ncol(SB)]
  uB <- xB[,-ncol(xA)]

  fit <- glm(cbind(yA, 1-yA)~xA, family = binomial)  # summary(fit)
  xi1 <- fit$coefficients[1]
  xi2 <- fit$coefficients[-1]
  theta <- nlminb(par0, ell)$par

  alpha <- theta[1]
  beta <- theta[-c(1, length(theta))]
  gamma <- theta[length(theta)]

  ### IPW
  piAA <- c(1/{1 + exp(alpha +  uA%*%beta + gamma*yA)})
  muIPW <- sum(yA/piAA)/sum(1/piAA)

  ### REG
  muB <- c(xi1 + xB%*%xi2)
  cB <- log((1 + exp(muB+gamma))/(1 + exp(muB)))
  piB <- 1/{1 + exp(alpha +  uB%*%beta + cB)}
  mB <- piB * plogis(muB) + (1 - piB)* plogis(muB+gamma)
  muREG <- sum(dB*mB)/sum(dB)

  ### AIPW
  muA <- c(xi1 + xA%*%xi2)
  cA <- log((1 + exp(muA+gamma))/(1 + exp(muA)))
  piA <- 1/{1 + exp(alpha + uA %*%beta + cA)}
  mA <- piA * plogis(muA) + (1 - piA)* plogis(muA+gamma)
  muAIPW <- sum((yA - mA)/piAA)/sum(1/piAA) + muREG

  list(theta = theta, xi = c(xi1, xi2), mu = c(muREG, muIPW, muAIPW))
}
