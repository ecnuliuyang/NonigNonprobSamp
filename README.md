# NonigNonprobSamp

The R package **NonigNonprobSamp** implements the pseudo-likelihood inference for a non-probability sample, utilizing auxiliary information from an existing reference probability sample. Additionally, three estimates are provided for the finite population mean, employing regression-based prediction, inverse probability weighting (IPW), and augmented IPW (AIPW) techniques.

+ MAR(): used to calculate the estimates under the ignorable participation mechanism.

- MNAR(): used to calculate the estimates under the nonignorable participation mechanism.

+ seMNAR(): used to calculate the standard errors of the estimators under the nonignorable participation mechanism.


# Usage

In the **R** software, the following codes are used to install the package:

install.packages("devtools")

library(devtools)

install_github("ecnuliuyang/NonigNonprobSamp")


# Case study 

To reproduce analyses of the ESPACOV survey data in Liu et al. (2024), we present the corresponding code as follows. 

+ Code to summarize the data sets
  
library(NonigNonprobSamp)

dim(sampA)

dim(sampB)


- Code to calculate the estimates of the population proportion of Spaniards experiencing good moods under the ignorable participation mechanism

mean(sampA[,11], na.rm = T) # Naive estimate

MAR (sampA, sampB)$mu # Chen et al. (2020)'s method


+ Code to calculate the estimates of the population proportion of Spaniards experiencing good moods and the standard errors under the nonignorable participation mechanism
  
estMNAR <- MNAR (sampA, sampB)

estMNAR$mu

se <- seNMAR (sampA, sampB, estMNAR)

se$seMu


- Code to calculate the estimates of the estimated regression coefficients in the nonignorable participation model

cbind(round(estMNAR$theta,3), round(se$seTheta,3), round(2*(1 - pnorm(abs(estMNAR$theta)/se$seTheta)),3))

+ Code to calculate the estimates of the estimated regression coefficients in the outcome regression model

cbind(round(estMNAR$xi,3), round(se$seXi,3), round(2*(1 - pnorm(abs(estMNAR$xi)/se$seXi)),3))


# Reference

Chen, Y., Li, P. and Wu, C. (2020). Doubly robust inference with nonprobability survey samples. *Journal of the American Statistical Association* **115**, 2011–2021.

Liu, Y., Yuan, M., Li, P. and Wu, C. (2024). Statistical inference with nonignorable non-probability survey samples. *arXiv*: 2410.02920.


For questions, comments or remarks about the code please contact Y. Liu at this email address <liuyangecnu@163.com>.
