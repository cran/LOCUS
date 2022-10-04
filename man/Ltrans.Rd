\name{Ltrans}
\alias{Ltrans}
\title{
Map a symmetric matrix into its upper triangle part. }
\description{
This function is to map the upper triganle part of a symmetric matrix into a vector. 
}
\usage{
Ltrans(X, d = TRUE)
}

\arguments{
  \item{X}{A symmetric matrix of dimentional V by V. }
  \item{d}{Whether to include the diagonal part of \code{X}.}
}

\value{
  A vector containing the upper triganle part of \code{X}. 
}

