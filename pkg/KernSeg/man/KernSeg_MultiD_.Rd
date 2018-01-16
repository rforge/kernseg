\name{KernSeg_MultiD_}
\alias{KernSeg_MultiD_}
\title{KernSeg MultiD }
\description{Core function to call C code for the mD case }
\usage{KernSeg_MultiD_(geno, Kmax, delta, min.size = 1, alpha = NULL, 
    kernel = "Gaussian", option = 0)}
\arguments{
  \item{geno}{a matrix n x m}
  \item{Kmax}{maximum number of changes}
  \item{delta}{parameter of the kernel (delta.cn, delta.baf)}
  \item{min.size}{minimum size of a segment -1 (1 = no constraint 2 = at least 2)}
  \item{alpha}{some other parameter}
  \item{kernel}{kernel ("L2", "Gaussian")}
  \item{option}{for Mattesson divergence 0 or 1 (as Kernseg or as Ecp divide by n or n-1 per segment)}
}



\author{Guillem Rigaill}





