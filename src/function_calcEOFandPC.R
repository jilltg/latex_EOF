## AUTHOR:    JILL THORR, M.S.
## NAME:      f.calcEOFandPC
## DATE:	    9 JAN 2023
## OBJECTIVE: CALCULATE THE EMPIRICAL ORTHOGONAL FUNCTION 
##			      (EOF) MODES (SPATIAL) AND PRINCIPAL COMPONENT 
##			      (PC) SCORES (TEMPORAL) FROM MATRIX DECOMPOSED
##		     	  VIA SINGLUAR VALUE DECOMPOSITION.
## OUTPUTS:
##	- EOF    = Normalized EOF modes (values between -1 to 1)
##	- PC     = Principal Component Amplitudes (unconstrained values)
##  - VarExp = Variance explained by the first 25 EOF modes 
##	   		   * Restricted to 25 modes to reduce computational time

f.calcEOFandPC <- function(svdresults, A) {

	# svdresults  = output from f.SVD 
	# A           = data matrix that was decomposed in SVD (e.g., detrendedResiduals)

	# Create dim vars
	ngrids     <- ncol(A) # number of grid cells decomposed in SVD analysis [e.g., nrow(LLID)]
	ntimesteps <- nrow(A) # number of time-steps decomposed in SVD analysis [e.g., nrow(d)]

	# Format SVD outputs
	A.SVD   <- svdresults$A.SDV 
	neof    <- svdresults$neof
	A.SVD$v <- as.matrix(svdresults$A.SVD$v)
	A.SVD$u <- as.matrix(svdresults$A.SVD$u)
	A.SVD$d <- as.matrix(svdresults$A.SVD$d)
	temp    <- A.SVD$d%*%t(A.SVD$v)
	  
	# Normalize by EOF using max absolute value for each EOF and scaling EOfs and PCs
	# 	Note: 	dim(EOF) = # grid cells x # eofs 
	#		       	dim(pc)  = # eofs x # time-steps
  EOF     <- matrix(nrow = ngrids, ncol = neof) 
	pc      <- matrix(nrow = neof,   ncol = ntimesteps)
	for(k in 1:neof){

		eof_max <- max(A.SVD$u[, k])
		eof_min <- min(A.SVD$u[, k])
	
		if(abs(eof_max) > abs(eof_min)){
		  SCALE  <- eof_max
		} else {
		  SCALE  <- eof_min
		}
	
		EOF[ ,k] <- A.SVD$u[, k]/SCALE 
		pc[k, ]  <- SCALE*temp[k, ]

	} # end for-loop

	
	# Calculate variance explained by each mode
	neof     <- dim(EOF)[2] 
	var_orig <- sd(as.matrix(A))^2                        # initial variance
	vareof   <- c(); B <- matrix(nrow = dim(A)[1], ncol = dim(A)[2])
	ntorun   <- ifelse(dim(EOF)[2] > 25, 25, dim(EOF)[2]) # only calculate for the first 25 modes
	for(k in 1:ntorun) {
	
	  for(j in 1:dim(EOF)[1]) {
	   	 B[,j]  <- EOF[j,k]*t(pc[k,])
	  } 
	
	  residx    <- A - B
	  var_resid <- sd(as.matrix(residx))^2
	  vareof[k] <- 100*(1-var_resid/var_orig)

	}

	return(list("EOF"   = EOF,
				"PC"    = pc,
				"VarExp"= vareof))

} # end function

	
