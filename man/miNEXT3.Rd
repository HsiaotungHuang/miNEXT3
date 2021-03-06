% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Mainfuns.r
\name{miNEXT3}
\alias{miNEXT3}
\title{Mixture diversity of 3 communities
\code{miNEXT} computes the mixture diversity of 3 communities.}
\usage{
miNEXT3(data, knots = 15, m1 = NULL, nboot = 0)
}
\arguments{
\item{data}{a matrix/data.frame.}

\item{knots}{an integer specifying the number of points of m1. Default to be 15.}

\item{m1}{a vector specifying the values of m1 where the mixture diversity will be computed.}

\item{nboot}{an integer specifying the number of bootstrap times to build confidence interval. 
Use 0 to skip bootstrap.}
}
\value{
a list of 2 components: \code{$Each} a table of diversity of single community; \code{$Mixture} a 
table of mixture diversity.
}
\description{
Mixture diversity of 3 communities
\code{miNEXT} computes the mixture diversity of 3 communities.
}
\examples{
\dontrun{
data(Abudata)
data1 <- Abudata$data
result <- miNEXT3(data1)
}
}
