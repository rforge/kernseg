
KernSeg_MultiD <- function(
### Function to run the dynamic programing with kernels (md case)
geno, 
### an n x m  matrix to be segmented with n the length of the signal
Kmax, 
### the maximum number of breaks
delta=1,
### parameter of the kernel (delta.cn, delta.baf)
min.size=1,
### minimum size of a segment -1 (1 = no constraint 2 = at least 2)
alpha = NULL, 
### some additionnal parameter for the kernel
kernel='Gaussian',
### kernel ("Gaussian", "L2")
option=0
### for KernSeg like Mattesson Divergence or Mattesson like Kernseg (n or n-1) 
){
	nRow <- nrow(geno)
	out        = KernSeg_MultiD_(geno, Kmax, alpha = alpha, delta, kernel=kernel, min.size=min.size, option=option)
	t.est      = t(sapply(1:Kmax,FUN=function(k){c(retour(out$M_ij, k), rep(0,Kmax-k))}))
	t.est[1,1] = nrow(geno)
	res        = list(J.est  = out$C_ij[, nRow] + sum(out$D_i),t.est  = t.est)
	return(res)
### return a list with J.est : a vector with the loss as a function of K) and t.est : a matrix with the position 
### of the changes 
}


KernSeg_CNBAF_ <- function(
### KernSeg for Cn and Baf core function
geno, 
### 2 column matrix (cn, baf)
Kmax, 
### maximum number of changes
delta,
### parameter of the kernel (delta.cn, delta.baf)
min.size=1
### minimum size of a segment -1 (1 = no constraint 2 = at least 2)
){
	nRow <- nrow(geno)
	d <- ncol(geno)
	D_i = double(nRow) 
      	### Appel C
	### Guassian
	A <- .C("KernDynProg_dD_CNBAF_Raccess", 
                x_i= as.double(t(geno)),
                n= as.integer(nRow),
                d= as.integer(d),
                K= as.integer(Kmax), 
                M_ij= integer(Kmax*nRow), 
                C_ij= double(Kmax*nRow), 
                delta= as.double(delta), 
		D_i = D_i,
		min_size = as.integer(min.size),
                PACKAGE="KernSeg.light")
      


	A$C_ij = matrix(A$C_ij, nc=nRow, byrow=TRUE)
	A$M_ij = matrix(A$M_ij, nc=nRow, byrow=TRUE)
	
	A$t.est      = t(sapply(1:Kmax,FUN=function(k){c(retour(A$M_ij, k), rep(0,Kmax-k))}))
	A$t.est[1,1] <- nrow(geno)
	return(A)
### output of the C code (internal use)
      }

KernSeg_CNBAF <- function(
### kernel segmentation for CN and BAF signal
geno, 
### a n x 2 matrix where the first column is cn and the second is baf 
Kmax, 
### maximum number of changes
delta,
### parameter of the kernel (delta.cn, delta.baf)
min.size=1
### minimum size of a segment -1 (1 = no constraint 2 = at least 2)
){
	nRow <- nrow(geno)
	out        = KernSeg_CNBAF_(geno, Kmax, delta, min.size=min.size)
	t.est      = t(sapply(1:Kmax,FUN=function(k){c(retour(out$M_ij, k), rep(0,Kmax-k))}))
	t.est[1,1] = nrow(geno)
	res        = list(J.est  = out$C_ij[, nRow] + sum(out$D_i),t.est  = t.est)

	return(res)
###
}


KernSeg_MultiD_ <- function(
### Core function to call C code for the mD case 
geno, 
### a matrix n x m
Kmax, 
### maximum number of changes
delta,
### parameter of the kernel (delta.cn, delta.baf)
min.size=1,
### minimum size of a segment -1 (1 = no constraint 2 = at least 2)
alpha=NULL, 
### some other parameter
kernel='Gaussian',
### kernel ("L2", "Gaussian")
option=0
### for Mattesson divergence 0 or 1 (as Kernseg or as Ecp divide by n or n-1 per segment)
){
	nRow <- nrow(geno)
	d <- ncol(geno)
      	### Appel C
        if(kernel == "L2"){
	A <- .C("KernDynProg_dD_L2_Raccess", 
			x_i= as.double(t(geno)),
			n= as.integer(nRow),
			d= as.integer(d), 
			K= as.integer(Kmax), 
			M_ij= integer(Kmax*nRow), 
			C_ij= double(Kmax*nRow), 
			delta= as.double(delta), 
			D_i= double(nRow), 
			min_size = as.integer(min.size),
			PACKAGE="KernSeg")
	}
	
	if(kernel == "Gaussian"){
	A <- .C("KernDynProg_dD_Gaussian_Raccess", 
			x_i= as.double(t(geno)),
			n= as.integer(nRow),
			d= as.integer(d), 
			K= as.integer(Kmax), 
			M_ij= integer(Kmax*nRow), 
			C_ij= double(Kmax*nRow), 
			delta= as.double(delta), 
			D_i= double(nRow),
			min_size = as.integer(min.size), 
			PACKAGE="KernSeg")
	}

	if(kernel == "DivMatt"){
	A <- .C("KernDynProg_dD_DivMatt_Raccess", 
			x_i= as.double(t(geno)),
			n= as.integer(nRow),
			d= as.integer(d), 
			K= as.integer(Kmax), 
			M_ij= integer(Kmax*nRow), 
			C_ij= double(Kmax*nRow), 
			delta= as.double(delta), 
			D_i= double(nRow), 
			min_size = as.integer(min.size),
			option=as.integer(option),
			PACKAGE="KernSeg")
	}
  
      
	A$C_ij = matrix(A$C_ij, nc=nRow, byrow=TRUE)
	A$M_ij = matrix(A$M_ij, nc=nRow, byrow=TRUE)
	
	return(A)
### output of the C code (internal use)
	} 

############

retour <- function(
### function to parse C code output and recover changes in i breaks
path,
### matrix Mij of C output 
i
### line to consider
){
	chaine <- integer(i)
	chaine[i] <- ncol(path)
	for(j in i:2)  chaine[j-1] <- path[j, chaine[j]]
	return(chaine)
### list of changes (internal use)
}

