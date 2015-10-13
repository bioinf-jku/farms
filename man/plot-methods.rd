\name{plot-methods}
\docType{methods}
\alias{plot}
\alias{plot-methods}
\alias{plot,INI_Calls,missing-method}
\title{Visualizes the distribution of informative and non-informatives genes}
\description{
  This function visualizes the distribution of informative and non-informative genes of a given instance of \code{\link[farms]{INI_Calls-class}}.
}
\section{Methods}{
\describe{
\item{\code{signature(x = "INI_Calls", y = "missing")}}{An instance of \code{\link[farms]{INI_Calls-class}}. }
}}
\usage{
          \S4method{plot}{INI_Calls,missing}(x)
          }
\arguments{
    	\item{x}{An instance of \code{\link[farms]{INI_Calls-class}}.}
          }
\value{\code{\link[Biobase]{exprSet-class}}}
\seealso{\code{\link{expFarms}}, \code{\link{qFarms}},\code{\link{lFarms}},\code{\link{INIcalls}},\code{\link{summary}}}
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

