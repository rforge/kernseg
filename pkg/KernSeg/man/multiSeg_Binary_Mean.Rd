\name{multiSeg_Binary_Mean}
\alias{multiSeg_Binary_Mean}
\title{multiSeg Binary Mean}
\description{Function to run binary segmentation ona m x n matrix}
\usage{multiSeg_Binary_Mean(geno, weights = NULL, Kmax)}
\arguments{
  \item{geno}{m x n matrix where m is the number of chanel and n the length of the data}
  \item{weights}{an optionnal weight vector}
  \item{Kmax}{the maximum number of changes}
}

\value{Output of the C code wit x_i the signal, weights, weights, K the max number of breaks
n the length of the signal, P the number of dimension, list of recovered change-points 
during the successive step of binary segmentation, J.est cost of the segmentations}

\author{Guillem Rigaill}





