\name{INI_Calls-class}
\title{Class INI_Calls}
\docType{class}
\alias{INI_Calls}
\alias{INI_Calls-class}
\description{This is a class representation for an \code{INI_calls-class} object. The \code{INI_calls-class} consists of two instances of \code{\link[Biobase]{exprSet-class}}, containing an informative \code{exprSet} and a non-informative \code{exprSet}.
}

\section{Objects from the Class}{
  Objects can be created using the function \code{\link{INIcalls}}.
}

\section{Slots}{
  \describe{
    \item{\code{I_Calls}:}{Object of class \code{"vector"} containing informative probe set names.}
    \item{\code{NI_Calls}:}{Object of class \code{"vector"}  containing non-informative probe set names.}
    \item{\code{I_Exprs}:}{Object of class  \code{exprSet-class} representing the informative \code{exprSet}.}
    \item{\code{NI_Exprs}:}{Object of class \code{exprSet-class} representing the non-informative \code{exprSet}.}
    \item{\code{varZX}:}{Object of class \code{"vector"} containing the INI-call value.}
    }
}

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
\author{Djork Clevert}
\keyword{classes}

