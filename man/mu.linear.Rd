% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mu.linear.R
\name{mu.linear}
\alias{mu.linear}
\title{Produce mean vector for linear}
\usage{
mu.linear(doses)
}
\arguments{
\item{doses}{a numeric vector of specified doses.}

\item{k}{a numeric value specifies the number of active doses.}
}
\value{
mean vector for linear
}
\description{
This function is a help function to automatrically produce mean vector for linear.
}
\examples{
k<-4
doses<-seq(0,k)
mu.linear(doses)

}

