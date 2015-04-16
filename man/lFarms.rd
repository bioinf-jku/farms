\name{lFarms}
\alias{lFarms}
\title{lFarms expression measure}
\description{
  This function converts an instance of \code{\link[affy:AffyBatch-class]{AffyBatch}}
  into an instance of \code{\link[Biobase]{exprSet-class}} using a factor analysis model 
  for which a Bayesian Maximum a Posteriori method optimizes the model parameters under 
  the assumption of Gaussian measurement noise. This function is a wrapper for \code{\link[affy]{expresso}} 
  and uses the function \code{\link[affy]{normalize.loess}} for array normalization.
}
\usage{
          lFarms(object, weight, mu, weighted.mean, laplacian, robust, correction, centering,spuriousCorrelation,...)
          }
\arguments{
    	\item{object}{An instance of \code{\link[affy:AffyBatch-class]{AffyBatch}}.}
	\item{weight}{Hyperparameter value in the range of [0,1]  which determines the influence of the prior. The default value is 0.5 } 
	\item{mu}{Hyperparameter value which allows to quantify different aspects of 
	potential prior knowledge. Values near zero assumes that most genes do not
	contain a signal, and introduces a bias for loading matrix elements near zero. Default value is 0}
 	\item{weighted.mean}{Boolean flag, that indicates whether a weighted mean or a least square fit is used to summarize the loading matrix. The default value is set to FALSE.}
 	\item{laplacian}{Boolean flag, indicates whether a Laplacian prior for the factor is employed or not. Default value is FALSE.}
 	\item{robust}{Boolean flag, that ensures non-constant results. Default value is TRUE.}
 	\item{correction}{Value that indicates whether the covariance matrix should be corrected for negative eigenvalues 
	which might emerge from the non-negative correlation constraints or not. Default = O  (means that no correction is done), 
	1 (minimal noise (0.0001) is added to the diagonal elements of the covariance matrix to force positive definiteness), 
	2 (Maximum Likelihood solution to compute the nearest positive definite matrix under the given non-negative correlation constraints of the covariance matrix)}
	\item{centering}{Indicates whether the data is "median" or "mean"  centered. Default value is "median".}
	 \item{spuriousCorrelation}{Numeric value in the range of [0,1] that quantifies the suppression of spurious correlation when using 
	 the Laplacian prior. Default value is 0 (no suppression). Note, that this parameter is only active when the laplacian parameter is set to TRUE.}
	\item{...}{other arguments to be passed to \code{\link[affy]{expresso}}.}
          }
\value{\code{\link[Biobase]{exprSet-class}}}
\details{
  	This function is a wrapper for \code{\link[affy]{expresso}}.}	
\seealso{\code{\link[affy]{expresso}}, \code{\link[farms]{expFarms}}, \code{\link[farms]{qFarms}}, \code{\link[affy]{normalize.loess}}}
\examples{
data(testAffyBatch)
eset <- qFarms(testAffyBatch)
}


\keyword{manip}

