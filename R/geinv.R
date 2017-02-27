# Returns the Moore-Penrose inverse of the argument"
geninv <- function(G){                             
  # Transpose if m < n"                               
  m <- nrow(G)
  n <- ncol(G)
  transpose <- FALSE               
   if (m<n){                                          
    transpose <- TRUE                                  
    A <- tcrossprod(G)                                          
    n <- m                                             
   }
   else {                                           
     A <- crossprod(G)
   }                                                  
   # Full rank Cholesky factorization of A            
   dA <- diag(A)
   tol <- min(dA[dA>0])*1e-9          
   L <- matrix(0, nrow=nrow(A), ncol=ncol(A))
   r <- 0                                             
   for (k in 1:n){                                    
     r <- r+1
     L[k:n,r] <- A[(k:n),k,drop=FALSE] - tcrossprod(L[(k:n),(1:(r-1))], L[k,(1:(r-1)), drop=FALSE])
     # Note: for r=1, the substracted vector is zero"    
     if (L[k,r]>tol){                                   
       L[k,r] <- sqrt(L[k,r])                             
       if (k<n){                                          
         L[(k+1):n,r] <- L[(k+1):n,r]/L[k,r]                
       }                                                  
     }
     else { 
       r <- r-1                                           
     }                                                  
   }
    L <- L[,1:r]                                      
    # Computation of the generalized inverse of G"      
    M <- solve(crossprod(L))
    if (transpose){                                    
       Y <- crossprod(G,L)%*%M%*%tcrossprod(M,L)                                   
    }
    else {                                           
     Y <- L%*%M%*%tcrossprod(tcrossprod(M,L), G)                                   
    }
    Y
}    
