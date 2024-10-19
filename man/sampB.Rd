\name{sampB}
\alias{sampB}
\docType{data}
\title{Probability sample from the ESPACOV survey data}

\description{
The survey was designed to assess the impact of the COVID-19 pandemic in Spain and was conducted from January 18 to 25, 2021, approximately one year into the pandemic. To gather data, a probability sample was obtained by sending Short Message Service invitations to mobile phone numbers generated at random. Missing data are omitted prior to our analyses. Refer to Rueda et al. (2023) for the original data set.
}

\usage{
  data(sampB)
}

\format{
Our study focuses primarily on investigating the self-assessment of mood among Spanish residents aged 18 and older during the COVID-19 crisis. The probability sample consists of 881 observations and 11 variables.

\describe{

\item{x11}{
A binary variable, equal to 1 if \emph{age} falls between 30 and 44, and 0 otherwise.}

\item{x12}{
A binary variable, equal to 1 if \emph{age} falls between 45 and 65, and 0 otherwise.}

\item{x13}{
A binary variable, equal to 1 if when \emph{age} exceeds 65 and 0 otherwise.}

\item{x2}{
A binary variable, equal to 1 when the \emph{education level} is either first or second degree, and 0 when it is higher than second degree.
}

\item{x3}{
A binary variable, equal to 1 for the female and 0 for the male.
}

\item{x4}{
A binary variable, equal to 1 if the score of government action is higher and 0 otherwise.
}

\item{x5}{
A binary variable, equal to 1 if the cost of obeying the policy is either nothing or a bit and 0 otherwise.
}

\item{x6}{
A binary variable, equal to 1 if he/she will get vaccinated next year or never and 0 otherwise.
}

\item{x7}{
A binary variable, equal to 1 if his/her self-assessment of social status is either good or very good and 0 otherwise.}

\item{x8}{
A binary variable, equal to 1 if his/her self-assessment of health is either good or very good and 0 otherwise.}

\item{weight}{
A numerical variable, representing the survey weights.}
}
}

\details{
This data is used as auxiliary information to infer the  probability of participating the nonprobability suvry sample.
}

\references{
Liu, Y., Yuan, M., Li, P. and Wu, C. (2024).
Statistical inference with nonignorable non-probability survey samples.
\emph{arXiv}: 2410.02920.
}

\source{
Rueda, M. D. M., Pasadas‐del‐Amo, S., Rodríguez, B. C., Castro‐Martín, L., & Ferri‐García, R. (2023). Enhancing estimation methods for integrating probability and nonprobability survey samples with machine‐learning techniques. An application to a Survey on the impact of the COVID‐19 pandemic in Spain. \emph{Biometrical Journal}, \bold{65}(2), 2200035.
}

\examples{
data(sampB)
apply(sampB, 2, mean)
}

\keyword{datasets}
