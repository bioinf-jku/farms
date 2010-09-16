\name{summary-methods}
\docType{methods}
\alias{summary}
\alias{summary-methods}
\alias{summary,INI_Calls-method}
\title{Summary of I/NI-calls}
\description{
  This function determinates the percentage of informative genes of a given instance of of \code{\link{INI_Calls-class}} 
  which has been summarized by \code{\link{expFarms}}, \code{\link{qFarms}} or \code{\link{lFarms}} before.}
\section{Methods}{
\describe{
\item{\code{signature(object = "INI_Calls")}}{An instance of \code{\link[farms]{INI_Calls-class}}. }
}}
\usage{\S4method{summary}{INI_Calls}(object,...)}
\arguments{
    	\item{object}{An instance of \code{\link{INI_Calls-class}} .}
\item{...}{extra arguments to pass to the respective function}
          }
\value{\code{\link[Biobase]{exprSet-class}}}
\seealso{\code{\link{expFarms}}, \code{\link{qFarms}},\code{\link{lFarms}},\code{\link{plot}},\code{\link{INIcalls}}}
\examples{
data(testAffyBatch)
eset <- expFarms(testAffyBatch, bgcorrect.method = "rma", pmcorrect.method = "pmonly", normalize.method = "constant") 
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

