\name{INIcalls-methods}
\docType{methods}
\alias{INIcalls}
\alias{INIcalls-methods}
\alias{INIcalls,ExpressionSet-method}
\title{Dimension reduction based on informative genes}
\description{
  This function generates an  instance of \code{\link[farms]{INI_Calls-class}} of given 
  which has been summarized by \code{\link{expFarms}}, \code{\link{qFarms}} or \code{\link{lFarms}} before, based on the informative genes.
}

\section{Methods}{
\describe{
\item{\code{signature(object = "ExpressionSet")}}{An instance of \code{\link[Biobase]{exprSet-class}}. }
}}
\usage{
\S4method{INIcalls}{ExpressionSet}(object)
          }
\arguments{
    	\item{object}{An instance of \code{\link[Biobase]{exprSet-class}}.}
          }
\value{\code{\link[Biobase]{exprSet-class}}}
\seealso{\code{\link{expFarms}}, \code{\link{qFarms}},\code{\link{lFarms}},\code{\link{INIcalls}}}
\examples{
data(testAffyBatch)
eset <- expFarms(testAffyBatch, bgcorrect.method = "none", pmcorrect.method = "pmonly", normalize.method = "constant") 
INIs <- INIcalls(eset)  # apply I/NI calls
summary(INIs)
plot(INIs) # draws a density plot of I/NI-calls
I_data <- getI_Eset(INIs) # affybatch containing only informative probe sets
NI_data <- getNI_Eset(INIs) # affybatch containing only non-informative probe sets
I_probes <- getI_ProbeSets(INIs) # vector containing only informative probe sets names
NI_probes <- getNI_ProbeSets(INIs) # vector containing only  non-informative probe sets names
}
\keyword{manip}
\keyword{methods}

