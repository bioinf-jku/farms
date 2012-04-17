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
		normalize.method = "quantiles", weight, mu,  weighted.mean, laplacian, robust, correction, centering=c("median","mean"), spuriousCorrelation, ...){
	
	if (missing(weight)){weight <- 0.5}
	
	if (missing(mu)){mu <- 0}
	
	if (missing(robust)){robust <- TRUE}
	
	if (missing(weighted.mean)){weighted.mean <- FALSE}
	
	if (missing(correction)){correction <- 0}
	
	if (missing(laplacian)){laplacian <- FALSE}
	
	if (missing(spuriousCorrelation)){spuriousCorrelation <- 0}
	
	centering <- match.arg(centering)
	
	res <- expresso(object, bgcorrect.method=bgcorrect.method, pmcorrect.method=pmcorrect.method, 
			normalize.method=normalize.method, summary.method = "farms", 
			summary.param=list(weight=weight, mu=mu,  weighted.mean=weighted.mean, robust=robust, correction=correction, laplacian=laplacian, centering=centering, spuriousCorrelation=spuriousCorrelation))
	
	return(res)
	
} 


qFarms<-function (object, weight, mu, weighted.mean, laplacian, robust, correction, centering=c("median","mean"), spuriousCorrelation,...){
	
	if (missing(weight)){weight <- 0.5}
	
	if (missing(mu)){mu <- 0}
	
	if (missing(robust)){robust <- TRUE}
	
	if (missing(weighted.mean)){weighted.mean <- FALSE}
	
	if (missing(correction)){correction <- 0}
	
	if (missing(laplacian)){laplacian <- FALSE}
	
	if (missing(spuriousCorrelation)){spuriousCorrelation <- 0}
	
	centering <- match.arg(centering)
	
	res <- expresso(object, bgcorrect.method = "none", pmcorrect.method = "pmonly", 
			normalize.method = "quantiles", summary.method = "farms", 
			summary.param=list(weight=weight, mu=mu, weighted.mean=weighted.mean, robust=robust, correction=correction, laplacian=laplacian, centering=centering, spuriousCorrelation=spuriousCorrelation))
	
	return(res)
	
}

lFarms<-function (object, weight, mu, weighted.mean, laplacian, robust, correction, centering=c("median","mean"), spuriousCorrelation, ...){
	
	if (missing(weight)){weight <- 0.5}
	
	if (missing(mu)){mu <- 0}
	
	if (missing(robust)){robust <- TRUE}
	
	if (missing(weighted.mean)){weighted.mean <- FALSE}
	
	if (missing(correction)){correction <- 0}
	
	if (missing(laplacian)){laplacian <- FALSE}
	
	if (missing(spuriousCorrelation)){spuriousCorrelation <- 0}
	
	centering <- match.arg(centering)
	
	res <- expresso(object, bgcorrect.method = "none", pmcorrect.method = "pmonly", 
			normalize.method = "loess", summary.method = "farms", 
			summary.param=list(weight=weight, mu=mu, weighted.mean=weighted.mean, robust=robust, correction=correction, laplacian=laplacian, centering=centering, spuriousCorrelation=spuriousCorrelation))
	
	return(res)
	
}


generateExprVal.method.farms <- function(probes, weight, mu,  cyc, tol, weighted.mean, robust, minNoise, correction, laplacian, centering=c("median","mean"),spuriousCorrelation, ...){
	
	if (missing(weight)){weight <- 0.5}	
	
	if (missing(mu)){mu <- 0}
	
	if (missing(tol)){tol <- 0.00001}
	
	if (missing(robust)){robust <- TRUE}
	
	if (missing(cyc)){cyc <- 30}
	
	if (missing(weighted.mean)){weighted.mean <- FALSE}
	
	if (missing(correction)){correction <- 0}
	
	if (missing(laplacian)){laplacian <- FALSE}
	
	if (missing(minNoise)){minNoise <- 0.0001}
	
	if (missing(spuriousCorrelation)){spuriousCorrelation <- 0}
	
	centering <- match.arg(centering)
	
	
	## probes - data matrix
	## weight - hyperparameter default (0.5)
	## mu - hyperparameter default (0)
	## scale - scaling parameter for quantiles- (1.5) and 
	## loess-normalization (2)
	## tol - termination tolerance (default = 0.00001)
	## cyc - maximum number of cycles of EM (default 100)
	## L - factor loadings
	## Psi - diagonal uniqueness noise matrix 
	
	a_old <- 1/0.5
	n_array <-  ncol(probes)
	n_probes <- nrow(probes)
	epsmin <- 1e-30
	
	
	
	if(n_array < 2){
		
		stop("Error: FARMS is a multi-array method and therefore not designed for single-array summarization!")
		
	}
	
	if(n_array < 4){
		
		message("Warning: FARMS is a multi-array method, therefore it is not recommended to apply FARMS for batch sizes smaller than 4 arrays!")
		
	}
	
	
	
	
	
	probes <- log2(probes)## log2-transform probe intensities
	
	
	if (centering=="median") {mean.probes <- rowMedians(probes)}
	
	if (centering=="mean") {mean.probes <- rowMeans(probes)}  ## calculate mean of probes
	
	
	centered.probes <- probes - mean.probes
	
	
	
	sd.probes <- sqrt(diag(crossprod(t(centered.probes))) / n_array) ## calculate sd of probes
	
	if(0 %in% sd.probes){
		
		index <- which(sd.probes < epsmin)
		
		sd.probes[index] <- 1	## avoiding division by zero
		
		probes <- probes / sd.probes ## standardize probes to variance 1
		
		x <- t(probes)
		
		if (centering=="median") {y_v <- rowMedians(probes)}
		
		if (centering=="mean") {y_v <- rowMeans(probes)}
		
		xmean <- matrix(y_v, n_array, n_probes, byrow = TRUE)
		
		X <- x - xmean  ## center data (0 mean)
		
		XX <- crossprod(X,X) / n_array
		
		diag(XX)[index] <- 1 ## avoiding division by zero
		
	}
	
	else{
		
		probes <- probes / sd.probes ## standardize probes to variance 1
		
		x <- t(probes)
		
		if (centering=="median") {y_v <- rowMedians(probes)}
		
		if (centering=="mean") {y_v <- rowMeans(probes)}
		
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
	
	Psi <- diagXX - L^2
	
	Psi[which(Psi < epsmin)] <- epsmin
	
	alpha <- weight * n_probes
	
	bbeta <- mu * alpha
	
	
	PsiOld <- Psi
	
	if(laplacian){    ## Variational approach if Laplacian prior was selected
		
		PsiL <- (1/Psi)*L
		
		a <- as.vector(1+crossprod(L,PsiL))
		
		bar <- PsiL/a
		
		beta <- t(bar)
		
		mu_ZX <- X%*%bar
		
		lapla <- 1/sqrt(mu_ZX^2)
		
		for (i in 1:cyc){
			
			## E Step
			
			PsiL <- (1/Psi)*L
			
			a <- 1/as.vector(as.vector(lapla)+crossprod(L,PsiL))
			
			mu_ZX <- X%*%PsiL*a
			
			EZZ <- mu_ZX^2+a
			
			## M Step
			
			sumXMU <- 1/n_array*crossprod(X,mu_ZX)
			
			L <- (sumXMU +Psi*bbeta)/(mean(EZZ)+Psi*alpha)
			
			L[which(L<0)] <- epsmin
			
			Psi <- diagXX-L*sumXMU+Psi*alpha*L*(mu-L)
			
			Psi[which(Psi < epsmin)] <- epsmin
			
			lapla <- 1/sqrt(mu_ZX^2 + a)
			
			if (spuriousCorrelation != 0) { 
				
				lapla[which(lapla < spuriousCorrelation)] <- spuriousCorrelation
				
			}
			
			if(max(abs(Psi - PsiOld)) / max(abs(PsiOld)) < tol) {
				
				break
				
			}
			
			PsiOld <- Psi
			
			#if (sqrt(sum((a_old - a)^2)) < tol){
			
			#	break
			
			#}
			
			#a_old <- a
			
		}
		
		c <- mu_ZX  ## hidden variable c - factor
		
		
		
		laplacian_SNR <- 1/as.vector(1+crossprod(L,PsiL)/as.vector(lapla))	
		
		
	}
	
	else{
		
		
		for (i in 1:cyc){
			
# E Step
			
			PsiL <- (1 / Psi) * L
			
			a <- as.vector(1 + crossprod(L, PsiL))
			
			bar <- PsiL / a
			
			beta <- t(bar)
			
			XXbeta <- XX %*% bar
			
			EZZ <- 1 - beta %*% L + beta %*% XXbeta
			
			t_XXbeta <- XXbeta + Psi * bbeta
			
			t_EZZ <- as.vector(EZZ) + Psi * alpha
			
			## M Step
			
			L <- t_XXbeta / t_EZZ
			
			Psi <- diagXX - XXbeta * L + Psi * alpha * L * (mu - L) 
			
			Psi[which(Psi < epsmin)] <- epsmin
			
			if(max(abs(Psi - PsiOld)) / max(abs(PsiOld)) < tol) {
				
				break
				
			}
			
			PsiOld <- Psi
			
			
			#if (sqrt((1/a_old - 1/a)^2) < tol){
			
			#	break
			
			#}
			
			#a_old <- a
			
		}
		
		c <- X %*% bar ## hidden variable c - factor
		
	}
	
	if (!laplacian){
		
		if(EZZ < epsmin){
			
			var_z_scale <- 1 ## avoiding division by zero
		}
		else{
			
			var_z_scale <- sqrt(EZZ)
			
		}
	}
	
	else{
		
		var_scale <- sd(as.vector(c))*(1-1/n_array)
		
		if(var_scale < epsmin){
			
			var_z_scale <- 1 ## avoiding division by zero
		}
		
		else{
			
			var_z_scale <- var_scale
			
		}
		
	}
	
	c <- c / as.vector(var_z_scale)
	
	L <- L * as.vector(var_z_scale)
	
	PsiL <- (1 / Psi) * L
	
	a <- as.vector(1 + crossprod(L,PsiL))
	
	SNR <- 1 / a ## INI-Call
	
	##	SIG <- as.vector(crossprod(L, diag(as.vector(1/Psi)))) %*% XX %*% diag(as.vector(1/Psi)) %*% L * a^-2 ## SIGNAL-Call
	
	signal_info <- numeric(length=n_array)
	
	
	if (weighted.mean){
		
		PsiLL <- ((1 / Psi) * L^2)
		
		sumPsiLL <- sum(PsiLL)
		
		propPsiLL <- PsiLL / sumPsiLL
		
		express <- as.vector(crossprod(L * sd.probes, propPsiLL)) * c + median(y_v * sd.probes)
		
	} 
	
	else
	
	{
		
		express <- median(L * sd.probes) * c + median(y_v * sd.probes)
		
	}
	
	
	
	if (laplacian){
		
		signal_info <- laplacian_SNR
		
	}
	
	else{
		
		signal_info[] <- SNR
		
	}
	
	if(robust && (sd(as.vector(express)) < 1e-5)){
		
		rowMedianX <- rowMedians(x)
		
		x <- x - rowMedianX
		
		x[which(x <= 0)] <- NA
		
		rowMeansX <- rowMeans(x,  na.rm = TRUE)
		
		rowMeansXCentered <-  rowMeansX - mean(rowMeansX)
		
		minVarExpress <- rowMeansXCentered / (sqrt(mean(rowMeansXCentered^2))+epsmin) * 1e-5
		
		express <- median(y_v * sd.probes) + minVarExpress
		
	}
	
	return(list(exprs=as.numeric(express),se.exprs=as.numeric(signal_info)))
}













