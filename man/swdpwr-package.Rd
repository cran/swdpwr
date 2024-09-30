\name{swdpwr-package}
\alias{swdpwr-package}
\alias{swdpwr}
\docType{package}
\title{
Power Calculation for Stepped Wedge Cluster Randomized Trials
}
\description{
This package includes a function \code{swdpower} that accounts for power calculation for stepped wedge cluster randomized trials.
}
\details{
\tabular{ll}{
Package: \tab swdpwr\cr
Type: \tab Package\cr
Version: \tab 1.9\cr
Date: \tab 2023-09-28\cr
License: \tab GPL (version 3)\cr
}
Previous literature and developement of software focused mainly on continuous outcomes and obtained approximation results for binary outcomes. This package implemented new methods of power calculation for stepped wedge designs with binary outcomes and also incorporated procedures for continuous outcomes. The function \code{swdpower} can accommodate both cross-sectional and cohort designs, binary and continuous outcomes, marginal (GEE) and conditional models (mixed effects model), three link functions (identity, log, logit links), with and without time effects (the default specification assumes no-time-effect) under exchangeable, nested exchangeable and block exchangeable correlation structures. Unequal numbers of clusters per sequence are also allowed. With this package, investigators can obtain more accurate calculation of statistical power, that will help a lot in the design and analysis of stepped wedge cluster randomized trials. The package also includes functions that perform calculations for the intra-cluster correlation coefficients based on the random effects variances as input variables for continuous and binary outcomes, respectively. Other supplementary documents can be found at: https://ysph.yale.edu/cmips/research/software/study-design-power-calculation/swdpwr/. The Shiny app for swdpwr can be accessed at: <https://jiachenchen322.shinyapps.io/swdpwr_shinyapp/>.

%~~ An overview of how to use the package, including the most important functions ~~
}
\author{
Jiachen Chen, Xin Zhou , Fan Li and Donna Spiegelman.

Maintainer:Jiachen Chen <jiachen.chen322@gmail.com>
}
\references{
Zhou X, Liao X, Kunz L M, et al. A maximum likelihood approach to power calculations for stepped wedge designs of binary outcomes[J]. Biostatistics, 2020, 21(1): 102-121.

Li F, Turner E L, Preisser J S. Sample size determination for GEE analyses of stepped wedge cluster randomized trials[J]. Biometrics, 2018, 74(4): 1450-1458.
}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{ stepped wedge }
\keyword{ cluster randomized trials }
\keyword{ sample size estimation }
\keyword{ cross-sectional designs}
\keyword{ cohort designs}
\keyword{ binary and continuous outcomes }
\keyword{ correlation structure }
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
%\examples{
%~~ simple examples of the most important functions ~~
%}
