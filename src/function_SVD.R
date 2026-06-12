## AUTHOR:    JILL THORR, M.S.
## NAME:      f.SVD
## DATE:	    9 JAN 2023
## OBJECTIVE: DECOMPOSE A DATA MATRIX VIA SINGULAR
##			      VALUE DECOPOSITION (SVD) 
## OUTPUTS:
##	- neof  = number of empirical orthogonal function modes
##	- A.SVD = decomposed components: $V, $U, $D
##				$V = 
##				$U = n_grids x n_time matrix 
##				$D = n_time  x n_time matrix

f.SVD   <- function(A) {
 
	# A = data matrix to be decomposed where rows = time, cols = space

	  neof      <- min(dim(A))  # neofs = the lesser of the two dims 	  
	  M         <- t(A)         # transpose so matrix is rows = space, cols = time
	  A.SVD     <- svd(M)
	  D         <- matrix(0, ncol = neof, nrow = neof)
	  diag(D)   <- (A.SVD$d)
	  A.SVD$d   <- D
	  
	  returnlist <- list("neof"   = neof,
					     "A.SVD"  = A.SVD)
	  
	  return(returnlist)

} #end function

