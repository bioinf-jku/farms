\name{getI_ProbeSets}
\docType{methods}
\alias{getI_ProbeSets}
\alias{getI_ProbeSets-methods}
\alias{getI_ProbeSets,INI_Calls-method}
\title{Method to generate a vector of informative probe set names}
\description{
  This function generates an  instance of   \code{vector-class}, that return a vector of  informative probe set names.
}
\section{Methods}{
\describe{
\item{\code{signature(object = "INI_Calls")}}{An instance of \code{\link[farms]{INI_Calls-class}}. }
}}
\usage{\S4method{getI_ProbeSets}{INI_Calls}(object)}
\arguments{
    	\item{object}{An instance of \code{\link[farms]{INI_Calls-class}}.}
          }
\value{\code{vector}}
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

