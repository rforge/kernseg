

multiSeg_Binary_Mean <- function(
### Function to run binary segmentation ona m x n matrix
geno, 
### m x n matrix where m is the number of chanel and n the length of the data
weights=NULL, 
### an optionnal weight vector
Kmax
### the maximum number of changes
){
	if(class(geno) == "matrix"){
	nRow <- nrow(geno)
	nCol <- ncol(geno)
	} else {
	nRow= length(geno)
	nCol=1
	}
	if(length(weights)!=nRow)  weights = rep(1, nRow)
	
            A <- .C("Call_BinSeg", 
		x_i= as.double((geno)),
		weights=as.double(weights),
		K= as.integer(Kmax),
		n= as.integer(nRow), 
		P= as.integer(nCol), 
		t.est= integer(Kmax),
		J.est = double(Kmax), 
		PACKAGE="KernSeg.light")
	    #A$Cost <- sum(geno^2) - sum(apply(geno, 2, sum)^2/nRow) + c(0, cumsum(A$RupturesCost))
	return(A)
### Output of the C code wit x_i the signal, weights, weights, K the max number of breaks
### n the length of the signal, P the number of dimension, list of recovered change-points 
### during the successive step of binary segmentation, J.est cost of the segmentations
}

###kernlab required for this
getApproxMat <-  function(
### function to get low rank approximation with kernlab
geno, 
### m x n matrix
referenceMat, 
### matrix de reference
KernFun
### kernel function as in kernlab
){
        kernMatI <- kernelMatrix(KernFun, referenceMat)
	kernMatA <- kernelMatrix(KernFun, geno, referenceMat)
	svdI<- svd(kernMatI)

	return(kernMatA %*% (svdI$u %*% diag((svdI$d)^-0.5)))
### return the approximate n x p matrix
}

cumsumSpec <- function(
### special function for restricted DP
geno
### m x n matrix
){
if(class(geno) == "matrix"){
	nRow <- nrow(geno)
	nCol <- ncol(geno)
	sumXs <- apply(geno, 2, cumsum)
	sumX2 <- c(0, cumsum(apply(geno^2, 1, sum)))
	} else {
	nRow= length(geno)
	nCol=1
	sumXs <- cumsum(geno)
	sumX2 <- c(0, cumsum(geno^2))
	}
A <- list()
A$sumXs = sumXs
A$sumX2 = sumX2
return(A)
### internal use
}



multiSeg_Binary_Kern <- function(
### binary segmentation for kernel using a low rank approx
geno, 
### m x n matrix
Kmax, 
### maximum number of changes
weights=NULL, 
### weights
KernFun=rbfdot(sigma = 0.1), 
### kernel function as in kernlab
referencePoint, 
### reference point for the approx
referenceMat=NULL
### reference mat
){

	if(class(geno) == "matrix"){
	nRow <- nrow(geno)
	nCol <- ncol(geno)
	} else {
	nRow= length(geno)
	nCol=1
	}

	### pre-calculation
	if(class(geno) == "matrix"){
	if(is.null(referenceMat))  referenceMat = geno[referencePoint, ]
	}
	
	if(class(geno) != "matrix"){
	if(is.null(referenceMat))  referenceMat = geno[referencePoint]
	}

	if(length(weights)!=nRow) weights=rep(1, nRow)

	ApproxMat <- getApproxMat(geno, referenceMat, KernFun)
	res <- multiSeg_Binary_Mean(ApproxMat, weights , Kmax)
	res$Approx <- ApproxMat

	cumSpec <- cumsumSpec(res$Approx)
        changes = c(sort(res$t.est), nrow(res$Approx))
        res_ <- multiSeg_RestrictedDynProg_Mean(cumSpec$sumXs, cumSpec$sumX2, c(0, cumsum(weights)), Kmax=Kmax, changes=changes)
	return(res_)
### return result after restricted DP on change points recovered after binary segmentation
}


multiSeg_Binary_Kern_withApprox <- function(
### binary segmentation of an m times n matrix followed by restricted DP
ApproxMat, 
### m x n matrix (m > 2)
Kmax, 
### maximum number of changes
weights=NULL
### weights
){

	if(is.null(weights)) weights=rep(1, nrow(ApproxMat))
	res <- multiSeg_Binary_Mean(ApproxMat, weights , Kmax)
	res$Approx <- ApproxMat

	cumSpec <- cumsumSpec(res$Approx)
        changes = c(sort(res$t.est), nrow(res$Approx))
        res_ <- multiSeg_RestrictedDynProg_Mean(cumSpec$sumXs, cumSpec$sumX2, c(0, cumsum(weights)), Kmax=Kmax, changes=changes)
	return(res_)
### return result after restricted DP on change points recovered after binary segmentation
}
