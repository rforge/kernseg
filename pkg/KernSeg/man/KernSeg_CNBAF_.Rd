\name{KernSeg_CNBAF_}
\alias{KernSeg_CNBAF_}
\title{KernSeg CNBAF }
\description{KernSeg for Cn and Baf core function}
\usage{KernSeg_CNBAF_(geno, Kmax, delta, min.size = 1)}
\arguments{
  \item{geno}{2 column matrix (cn, baf)}
  \item{Kmax}{maximum number of changes}
  \item{delta}{parameter of the kernel (delta.cn, delta.baf)}
  \item{min.size}{minimum size of a segment -1 (1 = no constraint 2 = at least 2)}
}

\value{output of the C code (internal use)}

\author{Guillem Rigaill}





