\name{multiSeg_Binary_Kern}
\alias{multiSeg_Binary_Kern}
\title{multiSeg Binary Kern}
\description{binary segmentation for kernel using a low rank approx}
\usage{multiSeg_Binary_Kern(geno, Kmax, weights = NULL, KernFun = rbfdot(sigma = 0.1), 
    referencePoint, referenceMat = NULL)}
\arguments{
  \item{geno}{m x n matrix}
  \item{Kmax}{maximum number of changes}
  \item{weights}{weights}
  \item{KernFun}{kernel function as in kernlab}
  \item{referencePoint}{reference point for the approx}
  \item{referenceMat}{reference mat}
}

\value{return result after restricted DP on change points recovered after binary segmentation}

\author{Guillem Rigaill}





