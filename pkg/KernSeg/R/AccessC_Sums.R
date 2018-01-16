
#### compute all block sums in gram matrix
#### from their we could get the cost matrix
SumsBlock_Gaussian <- function(geno, delta){
	nRow <- length(geno)
       
	A <- .C("sumsBlock_Gaussian", 
			x_i= as.double(geno),
			n= as.integer(nRow), 
			Sums= double(nRow*nRow),
			delta= as.double(delta), 
			PACKAGE="KernSeg.light")
	return(A)
}


