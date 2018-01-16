\name{multiSeg_RestrictedDynProg_Mean}
\alias{multiSeg_RestrictedDynProg_Mean}
\title{multiSeg RestrictedDynProg Mean}
\description{function to run restricted dp L2 (= DP on a restricted set of candidate change-points)}
\usage{multiSeg_RestrictedDynProg_Mean(sumXs, sumX2, sumWeights, 
    Kmax, changes)}
\arguments{
  \item{sumXs}{cumsum for X mxn}
  \item{sumX2}{cumsum of squares m x n}
  \item{sumWeights}{cumsum weights}
  \item{Kmax}{maximum number of changes}
  \item{changes}{position of changes}
}

\value{return result of C code with t.est and J.est (internal use)}

\author{Guillem Rigaill}





