\name{KernSeg_CNBAF}
\alias{KernSeg_CNBAF}
\title{KernSeg CNBAF}
\description{kernel segmentation for CN and BAF signal}
\usage{KernSeg_CNBAF(geno, Kmax, delta, min.size = 1)}
\arguments{
  \item{geno}{a n x 2 matrix where the first column is cn and the second is baf }
  \item{Kmax}{maximum number of changes}
  \item{delta}{parameter of the kernel (delta.cn, delta.baf)}
  \item{min.size}{minimum size of a segment -1 (1 = no constraint 2 = at least 2)}
}


\author{Guillem Rigaill}





