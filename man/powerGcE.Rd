\name{powerGcE}
\alias{powerGcE}
\title{powerGcE}
\description{
Evaluate gene by environment interaction with a binary outcome and a normally distributed environmental exposure.
}
\usage{
powerGcE(nCase = 407, nControl = 376, MAF = 0.49, meanE = 0, varE = 0.99, beta0 = -0.32, betaSNP = 0.17, betaE = 0.97, betaI = seq(-1, -0.75, by = 0.05), nSim = 1000, alpha = 5e-08, plot.output = TRUE, plot.name = "powerGcE.pdf", seed = 1)
}
\arguments{
  \item{nCase}{
is the number of cases
}
  \item{nControl}{
is the number of controls
}
  \item{MAF}{
is the minor allele frequency for the SNP
}
  \item{meanE}{
is the mean of the normally distributed environmental exposure

}
  \item{varE}{
is the variance of the normally distributed environmental exposure

}
  \item{beta0}{
For the binary outcome Y, the environmental exposure E, and the SNP, logit(P(Y=1))=Beta0+BetaSNP*SNP+BetaE*E+BetaI*SNP*E
}
  \item{betaSNP}{
For the binary outcome Y, the environmental exposure E, and the SNP, logit(P(Y=1))=Beta0+BetaSNP*SNP+BetaE*E+BetaI*SNP*E
}
  \item{betaE}{
For the binary outcome Y, the environmental exposure E, and the SNP, logit(P(Y=1))=Beta0+BetaSNP*SNP+BetaE*E+BetaI*SNP*E
}
  \item{betaI}{
For the binary outcome Y, the environmental exposure E, and the SNP, logit(P(Y=1))=Beta0+BetaSNP*SNP+BetaE*E+BetaI*SNP*E
}
  \item{nSim}{
is the number of simulations

}
  \item{alpha}{
is the alpha level, default=0.00000005
}
  \item{plot.output}{
if true, then a plot is outputted to the working directory.
}
  \item{plot.name}{
is the name of the plot.
}
  \item{seed}{
is set for reproducibility. 
}
}
\value{
The SNP is generated from a binomial distribution and the environmental exposure from a normal distribution. Then, the binary outcome is generated from a binomial distribution such that logit(P(Y=1))=Beta0+BetaSNP*SNP+BetaE*E+BetaI*SNP*E where E is the environmental exposure. We exmained the proportion of simulations where the p-value for the interaction term is less than alpha.
}
\author{
Sharon Lutz
}
\examples{
	powerGcE(nCase = 407, nControl = 376, MAF = 0.49, meanE = 0, varE = 0.99, beta0 = -0.32, betaSNP = 0.17, betaE = 0.97, betaI = seq(-1, -0.75, by = 0.05), nSim = 1000, alpha = 5e-08, plot.output = TRUE, plot.name = "powerGcE.pdf", seed = 1)
}
\keyword{powerGcE}
\keyword{power}
\keyword{GxE}
