##    Copyright (C) 2006 Djork-Arne Clevert (okko@clevert.de),
##       				 Sepp Hochreiter (hochreit@bioinf.jku.at),
##                       Klaus Obermayer (oby@cs.tu-berlin.de)
##    Berlin University of Technology,
##    Institute for Software Engineering and Theoretical Computer Science 
##    The software is maintained and developed by Djork-Arné Clevert. 
##    We offer a first implementation of the new 
##    ``Factor Analysis for Robust Microarray Summarization'' (FARMS) algorithm.
##    This program is free software; you can redistribute it and/or modify it under 
##    the terms of the GNU General Public License as published by the Free Software 
##    Foundation; either version 2 of the License, or (at your option) any later version. 
##    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
##    without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
##    See the GNU General Public License for more details.
##    If you use this library, please cite:
##
##    @article{SeppHochreiter02102006,
##		author = {Hochreiter, Sepp and Clevert, Djork-Arne and Obermayer, Klaus},
##		title = {{A new summarization method for Affymetrix probe level data}},
##		journal = {Bioinformatics},
##		volume = {},
##		number = {},
##		pages = {btl033},
##		doi = {10.1093/bioinformatics/btl033},
##		year = {2006},
##		URL = {http://bioinformatics.oxfordjournals.org/cgi/content/abstract/btl033v1},
##		eprint = {http://bioinformatics.oxfordjournals.org/cgi/reprint/btl033v1.pdf}
##		}

##	@article{INI-Calls:07,
##		author = {Willem Talloen and Djork-Arne Clevert and Sepp Hochreiter and Dhammika Amaratunga and Luc Bijnens and Stefan Kass and Hinrich W.H. Ghlmann},
##		title = {/NI-calls for the exclusion of non-informative genes: a highly effective filtering tool for microarray data},
##		journal = {Bioinformatics},
##		volume = {},
##		number = {},
##		pages = {btm478},
##		doi = {doi:10.1093/bioinformatics/btm478},
##		year = {2007},
##		URL = {http://bioinformatics.oxfordjournals.org/cgi/content/short/btm478v1},
##		eprint = {http://bioinformatics.oxfordjournals.org/cgi/reprint/btm478v1}
##		}

upDate.generateExprSet.methods(c(generateExprSet.methods(), "farms"))

upDate.express.summary.stat.methods(c(express.summary.stat.methods(), "farms"))





expFarms<-function(object, bgcorrect.method = "none", pmcorrect.method = "pmonly", 
normalize.method = "quantiles", weight, mu,  weighted.mean, laplacian, robust, correction,...){
	
	if (missing(weight)){weight <- 0.5}
	
	if (missing(mu)){mu <- 0}
	
	if (missing(robust)){robust <- TRUE}
	
	if (missing(weighted.mean)){weighted.mean <- TRUE}
	
	if (missing(correction)){correction <- 0}
	
	if (missing(laplacian)){laplacian <- FALSE}
	
    res <- expresso(object, bgcorrect.method=bgcorrect.method, pmcorrect.method=pmcorrect.method, 
					normalize.method=normalize.method, summary.method = "farms", 
					summary.param=list(weight=weight, mu=mu,  weighted.mean=weighted.mean, robust=robust, correction=correction, laplacian=laplacian))
	
    return(res)
	
} 


qFarms<-function (object, weight, mu, weighted.mean, laplacian, robust, correction, ...){
	
	if (missing(weight)){weight <- 0.5}
	
	if (missing(mu)){mu <- 0}
	
	if (missing(robust)){robust <- TRUE}
	
	if (missing(weighted.mean)){weighted.mean <- TRUE}
	
	if (missing(correction)){correction <- 0}
	
	if (missing(laplacian)){laplacian <- FALSE}
	
    res <- expresso(object, bgcorrect.method = "none", pmcorrect.method = "pmonly", 
					normalize.method = "quantiles", summary.method = "farms", 
					summary.param=list(weight=weight, mu=mu, weighted.mean=weighted.mean, robust=robust, correction=correction, laplacian=laplacian))
	
    return(res)
	
}

lFarms<-function (object, weight, mu, weighted.mean, laplacian, robust, correction,...){
	
	if (missing(weight)){weight <- 0.5}
	
	if (missing(mu)){mu <- 0}
	
	if (missing(robust)){robust <- TRUE}
	
	if (missing(weighted.mean)){weighted.mean <- TRUE}
	
	if (missing(correction)){correction <- 0}
	
	if (missing(laplacian)){laplacian <- FALSE}
	
    res <- expresso(object, bgcorrect.method = "none", pmcorrect.method = "pmonly", 
					normalize.method = "loess", summary.method = "farms", 
					summary.param=list(weight=weight, mu=mu, weighted.mean=weighted.mean, robust=robust, correction=correction, laplacian=laplacian))
	
    return(res)
	
}


generateExprVal.method.farms <- function(probes, weight, mu,  cyc, tol, weighted.mean, robust=FALSE, minNoise, correction, laplacian, ...){
	
	if (missing(weight)){weight <- 0.5}	
	
	if (missing(mu)){mu <- 0}
	
	if (missing(tol)){tol <- 0.00001}
	
	if (missing(robust)){robust <- TRUE}
	
	if (missing(cyc)){cyc <- 25}
	
	if (missing(weighted.mean)){weighted.mean <- TRUE}
	
	if (missing(correction)){correction <- 0}
	
	if (missing(laplacian)){laplacian <- FALSE}
	
	if (missing(minNoise)){minNoise <- 0.0001}
	
	
## probes - data matrix
## weight - hyperparameter default (0.5)
## mu - hyperparameter default (0)
## scale - scaling parameter for quantiles- (1.5) and 
## loess-normalization (2)
## tol - termination tolerance (default = 0.00001)
## cyc - maximum number of cycles of EM (default 100)
## L - factor loadings
## Ph - diagonal uniqueness matrix
	
	a_old <- 0.5
	n_array <-  ncol(probes)
	n_probes <- nrow(probes)
	

		
	if(n_array < 2){
			
		stop("Error: FARMS is a multi-array method and therefore not designed for single-array summarization!")
		
	}
		
	if(n_array < 4){
			
		message("Warning: FARMS is a multi-array method, therefore it is not recommended to apply FARMS for batch sizes smaller than 4 arrays!")
			
	}
		
		

	

	probes <- log2(probes)## log2-transform probe intensities
	
	mean.probes <- rowMeans(probes)  ## calculate mean of probes
	
	centered.probes <- probes - mean.probes
	
	sd.probes <- sqrt(diag(crossprod(t(centered.probes))) / n_array) ## calculate sd of probes
	
	if(0 %in% sd.probes){
		
		index <- which(sd.probes == 0)
		
		sd.probes[index] <- 1	## avoiding division by zero
		
		probes <- probes / sd.probes ## standardize probes to variance 1
		
		x <- t(probes)
		
		y_v <- colMeans(x)
		

		xmean <- matrix(y_v, n_array, n_probes, byrow = TRUE)
		
		X <- x - xmean  ## center data (0 mean)
		
		XX <- crossprod(X,X) / n_array
		
		diag(XX)[index] <- 1 ## avoiding division by zero
		
	}
	
	else{
		
		probes <- probes / sd.probes ## standardize probes to variance 1
		
		x <- t(probes)
		
		y_v <- colMeans(x)

		xmean <- matrix(y_v, n_array, n_probes, byrow = TRUE)
		
		X <- x - xmean  ## center data (0 mean)
		
		XX <- crossprod(X, X) / n_array
		
	}
	
	XX <- (XX + t(XX)) / 2 ## XX is now positive definit
	
	XX[which(XX < 0)] <- 0
	
	minEigenValues <- -1
	
	if(correction>=1){
		
		while(minEigenValues < 0){
			
    		eigen_XX <- eigen(XX)
			
			eigenValues_XX <- eigen_XX$values
			
			eigenVectors_XX <- eigen_XX$vectors
			
			minEigenValues <- min(eigenValues_XX)
			
    		if(correction<2){
				
				if(minEigenValues<minNoise){
					
					diag(XX)<-diag(XX)+(minNoise - minEigenValues)
					
				}
				
			}
			
    		else{
				
				if(minEigenValues<minNoise){
					
					eigenValues_XX[which(eigenValues_XX<minNoise)] <- minNoise
					
					XX <- eigenVectors_XX%*%diag(eigenValues_XX)%*%t(eigenVectors_XX)
					
				}
				
			}
			
		}
		
	}
	
	diagXX <- diag(XX)
	
	L <- sqrt(0.75 * diagXX) ## initialize factor loadings 
	
	Ph <- diagXX - L^2
	
    alpha <- weight * n_probes
	
	bbeta <- mu * alpha
	
	
	if(laplacian){    ## Variational approach if Laplacian prior was selected
		
		PsiL <- (1/Ph)*L
		
		a <- as.vector(1+crossprod(L,PsiL))
		
		bar <- PsiL/a
		
		beta <- t(bar)
		
		mu_ZX <- X%*%bar
		
		lapla <- 1/sqrt(mu_ZX^2)
		
		for (i in 1:cyc){
			
## E Step
			
			PsiL <- (1/Ph)*L
			
			L_robust <- L
			
			Ph_robust <- Ph
			
			a <- 1/as.vector(as.vector(lapla)+crossprod(L,PsiL))
			
			
			
			mu_ZX <- X%*%PsiL*a
			
			EZZ <- mu_ZX^2+a
			
## M Step
			
			sumXMU <- 1/n_array*crossprod(X,mu_ZX)
			
			L <- (sumXMU +Ph*bbeta)/(mean(EZZ)+Ph*alpha)
			
			L[which(L<0)] <- 0
			
			Ph <- diagXX-L*sumXMU+Ph*alpha*L*(mu-L)
			
			lapla <- 1/(mu_ZX^2 + a)^0.5
			
			if (sqrt(sum(a_old - a)^2) < tol){
				
				break
				
			}
			
			a_old <- a
			
		}
		
		c <- mu_ZX  ## hidden variable c - factor
		
		
		
		laplacian_SNR <- 1/as.vector(1+crossprod(L,PsiL)/as.vector(lapla))	
		
		
	}
	
	else{
		
		
		for (i in 1:cyc){
			
# E Step
			
			PsiL <- (1 / Ph) * L
			
			a <- as.vector(1 + crossprod(L, PsiL))
			
			L_robust <- L
			
			Ph_robust <- Ph
			
			if((1/a < 0.999) && robust){
				
				L_robust <- L
				
				Ph_robust <- Ph
				
			}
			
			bar <- PsiL / a
			
			beta <- t(bar)
			
			XXbeta <- XX %*% bar
			
			EZZ <- 1 - beta %*% L + beta %*% XXbeta
			
			t_XXbeta <- XXbeta + Ph * bbeta
			
			t_EZZ <- as.vector(EZZ) + Ph * alpha
			
## M Step
			
			L <- t_XXbeta / t_EZZ
			
			Ph <- diagXX - XXbeta * L + Ph * alpha * L * (mu - L) 
			
			if (sqrt(sum(1/a_old - 1/a)^2) < tol){
				
				break
				
			}
			
			a_old <- 1/a
			
		}
		
		
		
		c <- X %*% bar ## hidden variable c - factor
		
	}
	
	if (!laplacian){
		
		if(EZZ == 0){
			
			var_z_scale <- 1 ## avoiding division by zero
		}
		else{
			
			var_z_scale <- sqrt(EZZ)
			
		}
	}
	
	else{
		
		var_scale <- sd(c)*(1-1/n_array)
		
		if(var_scale == 0){
			
			var_z_scale <- 1 ## avoiding division by zero
		}
		
		else{
			
			var_z_scale <- var_scale
			
		}
		
	}
	
	c <- c / as.vector(var_z_scale)
	
	L <- L * as.vector(var_z_scale)
	
	PsiL <- (1 / Ph) * L
	
	a <- as.vector(1 + crossprod(L,PsiL))
	
	SNR <- 1 / a ## INI-Call
	
##	SIG <- as.vector(crossprod(L, diag(as.vector(1/Ph)))) %*% XX %*% diag(as.vector(1/Ph)) %*% L * a^-2 ## SIGNAL-Call
	
	
	signal_info <- numeric(length=n_array)
##	if (n_array >= 4){
##		signal_info[1] <- SNR
##		signal_info[2] <- SIG
##		signal_info[3] <- SIG * a^2
##		signal_info[4] <- i
##	}
##	if (n_array == 3){
##		signal_info[1] <- SNR
##		signal_info[2] <- SIG
##		signal_info[3] <- SIG * a^2
##	}
##	if (n_array == 2){
##		signal_info[1] <- SNR
##		signal_info[2] <- SIG
##	}
##	if (n_array == 1){
##		signal_info[1] <- SNR
##}
	
	
	if(robust && (SNR >= 0.999)){
		
		L <- L_robust
		
		Ph <- Ph_robust
		
		PsiL <- (1 / Ph) * L
		
		a <- as.vector(1 + crossprod(L, PsiL))
		
		bar <- PsiL / a
		
		beta <- t(bar)
		
		XXbeta <- XX %*% bar
		
		EZZ <- 1 - beta %*% L + beta %*% XXbeta
		
		c <- X %*% bar ## hidden variable c - factor
		
		if(EZZ == 0){
			
			var_z_scale <- 1 ## avoiding division by zero
			
		}
		
		else{
			
			var_z_scale <- sqrt(EZZ)
			
		}
		
		c <- c / as.vector(var_z_scale)
		
		L <- L * as.vector(var_z_scale)
		
	}
	
	
	
	if (weighted.mean){
		
		PsiLL <- ((1 / Ph) * L^2)
		
		sumPsiLL <- sum(PsiLL)
		
		propPsiLL <- PsiLL / sumPsiLL
		
		express <- as.vector(crossprod(L * sd.probes, propPsiLL)) * c + mean(y_v * sd.probes)
		
	} 
	
	else
	
	{
		
		express <- median(L * sd.probes) * c + mean(y_v * sd.probes)
		
	}
	
	
	
	if (laplacian){
		
		signal_info <- laplacian_SNR
		
	}
	
	else{
		
		signal_info[] <- SNR
		
	}
	
	return(list(exprs=as.numeric(express),se.exprs=as.numeric(signal_info)))
}














