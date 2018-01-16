\name{multiSeg_Binary_Kern_withApprox}
\alias{multiSeg_Binary_Kern_withApprox}
\title{multiSeg Binary Kern withApprox}
\description{binary segmentation of an m times n matrix followed by restricted DP}
\usage{multiSeg_Binary_Kern_withApprox(ApproxMat, Kmax, weights = NULL)}
\arguments{
  \item{ApproxMat}{m x n matrix (m > 2)}
  \item{Kmax}{maximum number of changes}
  \item{weights}{weights}
}

\value{return result after restricted DP on change points recovered after binary segmentation}

\author{Guillem Rigaill}





