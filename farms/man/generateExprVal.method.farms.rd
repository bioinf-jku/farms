\name{generateExprVal.method.farms}
\alias{generateExprVal.method.farms}
\title{Generate an expression value from the probes informations}
\description{Generate an expression from the probe}
\usage{
          generateExprVal.method.farms(probes, weight, mu,  cyc, tol, weighted.mean, robust, minNoise, correction, laplacian,  ...)
          }
\arguments{
   \item{probes}{a matrix of probe intesities with rows representing
    probes and columns representing
    samples. Usually \code{pm(probeset)} where \code{probeset} is a
    of class \code{\link[affy:ProbeSet-class]{ProbeSet}}}
	\item{weight}{Hyperparameter value in the range of [0,1]  which determines the influence of the prior. The default value is 0.5 } 
	\item{mu}{Hyperparameter value  which allows to quantify different aspects of 
	potential prior knowledge. A value near zero assumes that most genes do not
	contain a signal, and introduces a bias for loading matrix elements near zero. Default value is 0}
 	\item{cyc}{Value which determinates the maximum numbers of EM-Steps. Default value is set to number of arrays/2}
	\item{tol}{Value which determinates the termination tolerance. Convergence threshold is set to 1E-05.}
	\item{weighted.mean}{Boolean flag, that indicates wether a weighted mean or a least square fit is used to summarize the loading matrix. The default value is set to TRUE .}
 	\item{robust}{Boolean flag, that ensures non-constant results. Default value is TRUE.}
 	\item{minNoise}{Value, minimal noise assumption. Default value is 0.0001.}
 	 	\item{correction}{Value that indicates whether the covariance matrix should be corrected for negative eigenvalues 
	which might emerge from the non-negative correlation constraints or not. Default = O  (means that no correction is done), 
	1 (minimal noise (0.0001) is added to the diagonal elements of the covariance matrix to force positive definiteness), 
	2 (Maximum Likelihood solution to compute the nearest positive definite matrix under the given non-negative correlation constraints of the covariance matrix)}
	\item{laplacian}{Boolean flag, indicates whether a Laplacian prior for the factor is employed or not. Default value is FALSE.}
 		\item{...}{extra arguments to pass to the respective function}
    }

\value{
  A list containing entries:
  \item{exprs}{The expression values.}
  \item{se.exprs}{Estimate of the hidden variable.}
}
\seealso{
  \code{\link[affy]{generateExprSet-methods}},\code{\link[affy]{generateExprVal.method.playerout}},\code{\link[affy]{li.wong}}, \code{\link[affy]{medianpolish}}
}
\examples{

library(affy)
data(SpikeIn) ##SpikeIn is a ProbeSets
probes <- pm(SpikeIn)
exprs.farms <- generateExprVal.method.farms(probes)}

\keyword{manip}