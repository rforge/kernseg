multiSeg_RestrictedDynProg_Mean <- function(
### function to run restricted dp L2 (= DP on a restricted set of candidate change-points)
sumXs, 
### cumsum for X mxn
sumX2, 
### cumsum of squares m x n
sumWeights, 
### cumsum weights
Kmax, 
### maximum number of changes
changes
### position of changes
){
	
	    nRow =nrow(sumXs)
	    nCol = ncol(sumXs)
            A <- .C("Call_RestrictedDP", 
		sumXs_i= as.double((sumXs)),
		sumX2_i= as.double((sumX2)),
		sumW_i= as.double(sumWeights),
		n= as.integer(nRow), 
		P= as.integer(nCol), 
		changes = as.integer(changes-1),
		nChanges= length(changes),
		Kmax = as.integer(Kmax), 
		M_ij= integer(Kmax*length(changes)), 
		C_ij= double(Kmax*length(changes)), 
		PACKAGE="KernSeg.light")
		A$C_ij = matrix(A$C_ij, nrow=Kmax, byrow=T)
		A$M_ij = matrix(A$M_ij, nrow=Kmax, byrow=T)
		t.est      = t(sapply(1:Kmax,FUN=function(k){c(changes[retour(A$M_ij+1, k)], rep(0,Kmax-k))}))
		t.est[1,1] = length(A$sumX2_i)-1
		A$t.est = t.est 
		A$J.est = A$C_ij[, length(changes)]
	return(A)
### return result of C code with t.est and J.est (internal use)
}

