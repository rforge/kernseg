\name{KernSeg_MultiD}
\alias{KernSeg_MultiD}
\title{KernSeg MultiD}
\description{Function to run the dynamic programing with kernels (md case)}
\usage{KernSeg_MultiD(geno, Kmax, delta = 1, min.size = 1, 
    alpha = NULL, kernel = "Gaussian", option = 0)}
\arguments{
  \item{geno}{an n x m  matrix to be segmented with n the length of the signal}
  \item{Kmax}{the maximum number of breaks}
  \item{delta}{parameter of the kernel (delta.cn, delta.baf)}
  \item{min.size}{minimum size of a segment -1 (1 = no constraint 2 = at least 2)}
  \item{alpha}{some additionnal parameter for the kernel}
  \item{kernel}{kernel ("Gaussian", "L2")}
  \item{option}{for KernSeg like Mattesson Divergence or Mattesson like Kernseg (n or n-1) }
}

\value{return a list with J.est : a vector with the loss as a function of K) and t.est : a matrix with the position 
of the changes }

\author{Guillem Rigaill}





