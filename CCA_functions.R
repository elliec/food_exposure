options(scipen=999)

library(MASS)
library(tidyverse)
library(ggplot2)
library(GGally)
library(Matrix)
library(matrixcalc)
library(sf)

calcKmeans <- function(dat, fe, sa2s){
	set.seed(306411)
	clus_out <-dat
	vars <- NULL
	for (k in c(2:12)){
  		c <- kmeans(dat, centers=k, nstart=25,iter.max = 50)
  		clus_out <- cbind(clus_out, k=c$cluster)
  		vars <- rbind(vars, data.frame(totss = c$totss, withinss = c$tot.withinss, betweenss = c$betweenss))
	}

	clus_results <- data.frame(clus_out[,(ncol(clus_out)-10):ncol(clus_out)])
	names(clus_results) <- c(2:12)
	airclus = as.data.frame( clus_results ) %>%
  	pivot_longer(cols = 1:11,
               names_to = "max_clust",
               values_to = "cluster"
               )
	conf = xtabs( ~ cluster + as.integer(max_clust), airclus )

	names(clus_results) <- paste0("K",c(2:12))
	clus_results$SA2_MAIN16 <- fe$SA2_MAIN

	dat$SA2_MAIN16 <- fe$SA2_MAIN
	env.sp <- merge(sa2, dat, by="SA2_MAIN16")
	env.sp <- merge(env.sp, clus_results, by="SA2_MAIN16")
	
	vars$k <- c(2:12)
	n <- dim(dat)[1]
	vars$WV = vars$withinss/lead(vars$withinss)
	vars$CH = (vars$betweenss/(vars$k-1)) / (vars$withinss / (n-vars$k))

	return(list(env.sp=env.sp,conf=conf, vars=vars))
}

#function to calculate sqrtcov1inv  and sqrtcov2inv  
sinv = function(X1,X2){
	Xall = cbind(X1, X2)

	d1 = dim(X1)[2]
	d2 = dim(X2)[2]
	d1p = d1+1
	d = d1+d2

	# calculate covs
	covall = cov(Xall)
	cov1 = covall[1:d1,1:d1]
	r1 = rankMatrix(cov1)  # you need to use library(Matrix) for this
	cov2 = covall[d1p:d,d1p:d]
	r2 = rankMatrix(cov2)
	cov12 = covall[1:d1,d1p:d]
	
	if(is.singular.matrix(cov1)) {
		# do eigen and all that
		spectral1 = eigen( cov1 )
   		V1 = spectral1$vectors
   		lambda1 = spectral1$values
 		sqrtcov1inv = V1[,1:r1] %*% diag( sqrt( 1/lambda1[1:r1] ) )%*% t( V1[,1:r1]  ) 
		cov1inv  = V1[,1:r1] %*% diag( lambda1[1:r1]^(1/2)) %*% t( V1[,1:r1] ) 
	} else {
		spectral1 = eigen( cov1 )
		V1 = spectral1$vectors
		lambda1 = spectral1$values
		sqrtcov1inv  = V1 %*% diag( sqrt( 1/lambda1 ) ) %*% t( V1 )
		cov1inv  = V1 %*% diag( lambda1^(1/2)) %*% t( V1 )

	}

	if(is.singular.matrix(cov2)) {
		# do eigen and all that
		spectral2 = eigen( cov2 )
   		V2 = spectral2$vectors
   		lambda2 = spectral2$values
		sqrtcov2inv = V2[,1:r2] %*% diag( sqrt( 1/lambda2[1:r2] ) )%*% t( V2[,1:r2]  )
		cov2inv  = V2 %*% diag( lambda2^(1/2) ) %*% t( V2 )

	} else {
		spectral2 = eigen(cov2)
		V2 = spectral2$vectors
		lambda2 = spectral2$values
		sqrtcov2inv = V2 %*% diag( sqrt( 1/lambda2) ) %*% t( V2 )
		cov2inv  = V2 %*% diag( lambda2^(1/2) ) %*% t( V2 )
	}
	
	return(list(sqrtcov1inv, sqrtcov2inv, cov12, cov1inv,cov2inv))
}

#function to calculate norms of each column
colnorm <- function(x) {
	y <- x^2
	y <- sum(y)
	y <- x/sqrt(y)
	return(y)
}

calcCCA <- function(X1, X2) {
	n1 <- names(X1)
	n2 <- names(X2)

	X1 <- scale(X1)
	X2 <- scale(X2)

	#invert matrices
	s12inv <- sinv(X1,X2)

	#calculate approximate Chat
	AChat <- s12inv[[1]] %*% s12inv[[3]] %*% s12inv[[2]]

	#use svd to get correlations, calculate scores and loadings
	svdc <- svd( AChat) #correlations
	sings <- svdc$d

	#calculate weights phi
	phi <- t(svdc$u) %*% s12inv[[1]]
	phi <- t(phi)
	phin <- apply(phi,2, FUN=colnorm)
	row.names(phin) <- n1

	#calculate correlation coefficients from Inge's book
	x1u <- data.frame(s12inv[[4]] %*% phi)
	row.names(x1u) <- n1
	names(x1u) <- paste0("U", c(1:dim(X2)[2]))
	x2u <- data.frame(s12inv[[2]] %*% t(s12inv[[3]]) %*% phi)
	row.names(x2u) <- n2
	names(x2u) <- paste0("U", c(1:dim(X2)[2]))

	#calculate weights psi
	psi <- t(svdc$v) %*% s12inv[[2]]
	psi <- t(psi)
	psin <- apply(psi,2, FUN=colnorm)
	row.names(psin) <- n2

	#calculate correlation coefficients from Inge's book
	x1v <- data.frame(s12inv[[1]] %*% s12inv[[3]] %*% psi)
	row.names(x1v) <- n1
	names(x1v) <- paste0("V", c(1:dim(X2)[2]))

	x2v <- data.frame(s12inv[[5]] %*% psi)
	row.names(x2v) <- n2
	names(x2v) <- paste0("V", c(1:dim(X2)[2]))

	## calculate U scores
	U <- t(svdc$u) %*% s12inv[[1]] %*% t(scale(X1,scale=F))
	U <- t(U)

	#calculate V scores
	V <- t(svdc$v) %*% s12inv[[2]] %*% t(scale(X2,scale=F))
	V <- t(V)

	#loadings
	cor.u <- cor(X1,U)
	cor.v <- cor(X2,V)

	#calculate cross-loadings
	cor.x1v <- cor(X1,V)
	cor.x2u <- cor(X2,U)

	return(list(cor=sings, U=U, phi=phin, V=V, psi=psin, cor.u=cor.u, cor.v=cor.v, cor.x1v = cor.x1v, cor.x2u = cor.x2u,
		x1u=x1u,x2u=x2u,x1v=x1v,x2v=x2v))
}

## function test significance

Tk = function( k, n, d1, d2, vv ){
  	# vv is vector of singular values of Chat (returned by cancor in the component $cor)
	rr = length( vv )  
  	Tkout = - ( n - ( d1 + d2 + 3 )/2 ) * log( prod ( 1 - vv[ (k + 1) : rr ]^2 ) )
  	# compare with chisq on (d1 - k ) * (d2 - k ) dof
  	dof = ( d1 - k ) * ( d2 - k )
  	pval = pchisq( Tkout, df = dof, lower.tail = FALSE )
  	list( Tkout = Tkout, pval = pval )
}
