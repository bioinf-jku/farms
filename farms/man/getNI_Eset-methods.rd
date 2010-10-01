\name{getNI_Eset-methods}
\docType{methods}
\alias{getNI_Eset}
\alias{getNI_Eset-methods}
\alias{getNI_Eset,INI_Calls-method}


\title{Method to generate an ExpressionSet of  non-informative genes}
\description{
  This function generates an  instance of  \code{\link[Biobase]{exprSet-class}}, that contains only  non-informative probe sets.
}
\section{Methods}{
\describe{
\item{\code{signature(object = "INI_Calls")}}{An instance of \code{\link[farms]{INI_Calls-class}}. }
}}
\usage{\S4method{getNI_Eset}{INI_Calls}(object)}
\arguments{
    	\item{object}{An instance of \code{\link[farms]{INI_Calls-class}}.}
          }
\value{\code{\link[Biobase]{exprSet-class}}}
\seealso{\code{\link{expFarms}}, \code{\link{qFarms}},\code{\link{lFarms}},\code{\link{INIcalls}},\code{\link{summary}}}
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

