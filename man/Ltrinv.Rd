\name{Ltrinv}
\alias{Ltrinv}
\title{
Inverse function of \code{Ltrans} to map back a vector into a symmetric matrix. }
\description{
This function is the inverse function of \code{Ltrans}, which is to map a vector back to a symmetric matrix. 
}
\usage{
Ltrinv(x, V, d = TRUE)
}

\arguments{
  \item{x}{A vector to convert to a matrix, which is of length p. }
  \item{V}{Dimension of the matrix which \code{x} is converted to. }
  \item{d}{Whether diagonal is kept in \code{x} or not. }
}

\value{
  A symmetric matrix whose upper triangle part is \code{x}. 
}

