##' Doubly robust inference with nonprobability survey samples
##'
##' @description Function to apply the pseudo-likelihood method in the assumption of a nonignorable participation mechanism.
##'
##' @param SA A matrix, the non-probability  sample. The last column indicates the outcome variable and the other columns represent the covariates. Additionally, the second-to-last column serves as the instrumental variable.
##' @param SB A matrix, the reference probability sample. The last column indicates the survey weights and the other columns represent the covariates. Additionally, the second-to-last column serves as the instrumental variable.
##' @param estMNAR A list, the output object of the function *NMAR()*.
##' @param secOrderProb A character. "Hajek" indicates that the second-order selection probability is estimated using the Hajek (1964)'s approximation. Otherwise, it is estimated through simple random sampling without replacement.
##' @param N0 A numeric, the size of the finite population. If *NULL*, it is set to the sum of weights from the reference probability sample.
##'
##' @return A list with three elements:
##' \itemize{
##'   \item seTheta, the standard errors of the maximum pseudo-likelihood estimates of the parameters in the participation probability model.
##'   \item seXi, the standard errors of the maximum likelihood estimates of the parameters in the outcome regression model.
##'   \item seMu, the standard errors of the regression, IPW, and AIPW estimates of the population mean.
##' }
##'
##' @references
##' Hajek, J. (1964). Asymptotic theory of rejective sampling with varying probabilities from a finite population. \emph{Annals of Mathematical Statistics} **35**, 1491–1523.
##'
##' Liu, Y., Yuan, M., Li, P. and Wu, C. (2024).
##' Statistical inference with nonignorable non-probability survey samples.
##' \emph{arXiv}: 2410.02920.
##'
##' @importFrom stats plogis
##' @importFrom samplingVarEst Pkl.Hajek.s
##'
##' @export
##'
### Calculate the standard errors of our proposed estimates
seMNAR <- function (SA, SB, estMNAR, secOrderProb = "Hajek", N0=NULL) {
  muREG <- estMNAR$mu[1]
  muIPW <- estMNAR$mu[2]
  muAIPW <- estMNAR$mu[3]
  theta <- estMNAR$theta
  alpha <- theta[1]
  beta <- theta[-c(1, length(theta))]
  gamma <- theta[length(theta)]

  xi1 <- estMNAR$xi[1]
  xi2 <- estMNAR$xi[-1]

  nA <- nrow(SA)
  nB <- nrow(SB)

  yA <- c(SA[,ncol(SA)])
  xA <- SA[, -ncol(SA)]
  uA <- xA[, -ncol(xA)]

  dB <- SB[,ncol(SB)]
  xB <- SB[,-ncol(SB)]
  uB <- xB[, -ncol(xB)]

  muA <- c(xi1 + xA%*%xi2)
  muB <- c(xi1 + xB%*%xi2)

  qA <- c(uA%*%beta)
  qB <- c(uB%*%beta)
  piAA <- c(1/{1 + exp(alpha +  qA + gamma*yA)})
  NA_hat <- sum(1/piAA)

  cA1 <- plogis(muA)
  cA0 <- plogis(muA + gamma)
  cB1 <- plogis(muB)
  cB0 <- plogis(muB + gamma)
  cA <- log((1 - cA1)/(1 - cA0))
  cB <- log((1 - cB1)/(1 - cB0))

  piA <- c(1/( 1 + exp(alpha + qA + cA) ))
  piB <- c(1/( 1 + exp(alpha + qB + cB) ))
  mA <- piA *cA1 + (1 - piA)*cA0
  mB <- piB *cB1 + (1 - piB)*cB0

  hiAA <- cbind(1, uA, yA)
  hiA <- cbind(1, uA, cA0)
  hiB <- cbind(1, uB, cB0)

  ### variance of the IPW estimator
  V12 <- c(t(hiAA)%*%((1/piAA - 1)*(yA - muIPW))/NA_hat)
  V22 <- - t(hiA)%*%diag(1 - piA)%*%hiA/sum(1/piA)
  V23 <- - t(hiA) %*% diag((1 - piA)*(cA0-cA1)) %*% cbind(1, xA)/sum(1/piA)
  V33 <- t(cbind(1, xA)) %*% diag(- cA1 *(1-cA1)) %*% cbind(1, xA)/NA_hat

  V12e <- t(t(hiA)%*%((1 - piA)*(cA0 - cA1))/sum(1/piA) +
              c(rep(0,length(theta)-1),sum((1/piA - 1)*cA0*(1-cA0))/sum(1/piA)))


  hN <- sum((yA - mA)/piAA)/NA_hat
  V12a <- c(t(hiAA)%*%((1/piAA - 1)*(yA - mA - hN))/NA_hat)
  vb_theta <- hiB %*% solve(V22) * piB * dB
  vb_ipw <- hiB%*%t(V12%*%solve(V22))*piB*dB
  vb_reg <- (mB - muREG - hiB%*%t(V12e%*%solve(V22))*piB)*dB
  vb_aipw <- (mB - sum(mA/piA)/sum(1/piA) - piB*(hiB%*%solve(V22)%*%V12a))*dB
  if(is.null(N0)) N0 <- sum(dB)

  if (secOrderProb== "Hajek") {
    ### used in María del Mar Rueda (2022, Biometrical Journal)
    pijB <- Pkl.Hajek.s(1/dB)
    VB_theta <- t(vb_theta) %*% (1 - as.matrix(1/dB) %*%
                                   t(as.matrix(1/dB))/pijB) %*% vb_theta/N0
    VB_ipw <- c(t(vb_ipw) %*% (1 - as.matrix(1/dB) %*%
                                 t(as.matrix(1/dB))/pijB) %*% vb_ipw/N0)
    VB_reg <- t(vb_reg) %*% (1 - as.matrix(1/dB) %*%
                               t(as.matrix(1/dB))/pijB) %*% vb_reg/N0

    VB_aipw <- t(vb_aipw) %*% (1 - as.matrix(1/dB)%*%
                                 t(as.matrix(1/dB))/pijB)%*%vb_aipw/N0
  }else{
    VB_theta <- t(vb_theta) %*% (diag(nB*(N0-1)/(nB-1)/N0-nB/N0, nrow(hiB)) +
                                   matrix(1- nB*(N0-1)/(nB-1)/N0,nrow(hiB),
                                          nrow(hiB))) %*% vb_theta/N0
    VB_ipw <- c(t(vb_ipw) %*% (diag(nB*(N0-1)/(nB-1)/N0-nB/N0, nrow(hiB)) +
                                 matrix(1- nB*(N0-1)/(nB-1)/N0,nrow(hiB),
                                        nrow(hiB))) %*% vb_ipw/N0)
    VB_reg <- t(vb_reg) %*% (diag(nB*(N0-1)/(nB-1)/N0-nB/N0, nrow(hiB)) +
                               matrix(1- nB*(N0-1)/(nB-1)/N0,nrow(hiB),
                                      nrow(hiB))) %*% vb_reg/N0
    VB_aipw <- t(vb_aipw) %*% (diag(nB*(N0-1)/(nB-1)/N0-nB/N0, nrow(hiB)) +
                                 matrix(1- nB*(N0-1)/(nB-1)/N0,nrow(hiB),
                                        nrow(hiB)))%*%vb_aipw/N0
  }

  seXi <- solve(V33) %*% t(cbind(1, xA)) %*% diag((1-piAA)*(yA-cA1)^2) %*%
    cbind(1, xA) %*% solve(V33)/NA_hat
  seXi <- sqrt(diag(seXi/N0))

  ind <- hiA %*% solve(V22) + diag(yA-cA1) %*% cbind(1, xA) %*%
    t(solve(V22)%*%V23%*%solve(V33))
  seTheta <- t(ind) %*% diag(1 - piAA) %*% ind /NA_hat + VB_theta
  seTheta <- sqrt(diag(seTheta/N0))

  seIPW <- sqrt((sum(((yA - muIPW)/piAA + hiA %*% t(V12%*%solve(V22)) +
                        diag(yA-cA1) %*% cbind(1, xA) %*% t(V12%*%solve(V22)%*%V23%*%solve(V33)))^2*
                       (1 - piAA))/NA_hat + VB_ipw)/N0)

  ### variance of the regression estimator
  V13e <- t(t(cbind(1, xA))%*%( (1 - piA)*(cA0 - cA1)^2 +
                                  cA1*(1-cA1) + cA0*(1-cA0) * (1/piA -1)
  )/sum(1/piA))

  seREG <- sqrt((sum((hiA %*% t(V12e%*%solve(V22)) +
                        diag(yA-cA1) %*% cbind(1, xA) %*%
                        t((V12e%*%solve(V22)%*%V23 -V13e)%*%solve(V33)))^2*
                       (1 - piAA))/NA_hat + VB_reg)/N0)

  ### variance of the AIPW estimator
  varAIPW <- sum(( (yA - mA - hN)/piAA + hiA%*%solve(V22)%*%V12a +
                     (cbind(1,xA)*(yA-cA1))%*%
                     solve(V33)%*%t(V23)%*%
                     solve(V22)%*%V12a )^2*(1-piAA))/NA_hat + VB_aipw
  seAIPW <- sqrt(varAIPW/N0)

  list(seMu = c(seREG, seIPW, seAIPW), seXi = seXi, seTheta = seTheta)
}
