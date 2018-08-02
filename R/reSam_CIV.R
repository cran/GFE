#' Gross flows variance estimation.
#'
#' @description Gross flows variance estimation according to resampling method (Bootstrap or Jackknife).
#' @param sampleBase An object of class data.frame or data.table containing the sample selected to estimate the gross flows.
#' @param nRepBoot The number of replicates for the bootstrap method.
#' @param model A character indicating the model that will be used for estime the gross flows. The available models are: 'I','II','III','IV'.
#' @param niter The number of iterations for the \eqn{\eta_{i}} and \eqn{p_{ij}} model parameters.
#' @param type A character indicating the resampling method \emph{("Bootstrap" or "Jackknife")}
#' @param colWeights The data colum name containing the sampling weights to be used on the fitting process.
#' @param nonrft a logical value indicating the non response for the first time. 
#' @return \code{reSamGF} returns a list that contains the variance of each parameter of the selected model.
#' @details
#' The resampling methods for variance estimation are:
#' \describe{
#'      \item{Bootstrap:}{ This technique allows to estimate the sampling distribution of almost any statistic by using random sampling methods. Bootstrapping is the practice of estimating properties of an statistic (such as its variance) by measuring those properties from it's approximated sample.}
#'      \item{Jackknife:}{ The jackknife estimate of a parameter is found by systematically leaving out each observation from a dataset and calculating the estimate and then finding the average of these calculations. Given a sample of size \emph{n}, the jackknife estimate is found by aggregating the estimates of each \emph{n-1}-sized sub-sample.}
#' }
#' @examples
#' library(TeachingSampling)
#' library(data.table)
#' # Colombia's electoral candidates in 2014
#' candidates_t0 <- c("Clara","Enrique","Santos","Martha","Zuluaga","Blanco", "NoVoto")
#' candidates_t1 <- c("Santos","Zuluaga","Blanco", "NoVoto")
#' 
#' N 	   <- 100000
#' nCanT0  <- length(candidates_t0)
#' nCanT1  <- length(candidates_t1)
#' 
#' # Initial probabilities 
#' eta <- matrix(c(0.10, 0.10, 0.20, 0.17, 0.28, 0.1, 0.05),
#' 			   byrow = TRUE, nrow = nCanT0)
#' # Transition probabilities
#' P <- matrix(c(0.10, 0.60, 0.15, 0.15,
#' 			  0.30, 0.10, 0.25,	0.35,
#' 			  0.34, 0.25, 0.16, 0.25,
#' 			  0.25,	0.05, 0.35,	0.35,
#' 			  0.10, 0.25, 0.45,	0.20,
#' 			  0.12, 0.36, 0.22, 0.30,
#' 			  0.10,	0.15, 0.30,	0.45),
#' 	 byrow = TRUE, nrow = nCanT0)
#' 
#' citaMod <- matrix(, ncol = nCanT1, nrow = nCanT0)
#' row.names(citaMod) <- candidates_t0
#' colnames(citaMod) <- candidates_t1
#' 
#' for(ii in 1:nCanT0){ 
#'  citaMod[ii,] <- c(rmultinom(1, size = N * eta[ii,], prob = P[ii,]))
#' }
#' 
#' # # Model I
#' psiI   <- 0.9
#' rhoRRI <- 0.9
#' rhoMMI <- 0.5
#' 
#' citaModI <- matrix(nrow = nCanT0 + 1, ncol = nCanT1 + 1)
#' rownames(citaModI) <- c(candidates_t0, "Non_Resp")
#' colnames(citaModI) <- c(candidates_t1, "Non_Resp")
#' 
#' citaModI[1:nCanT0, 1:nCanT1] 		 <- P * c(eta) * rhoRRI * psiI  
#' citaModI[(nCanT0 + 1), (nCanT1 + 1)]  <- rhoMMI * (1-psiI) 
#' citaModI[1:nCanT0, (nCanT1 + 1)]  	 <-  (1-rhoRRI) * psiI * rowSums(P * c(eta))
#' citaModI[(nCanT0 + 1), 1:nCanT1 ] 	 <-  (1-rhoMMI) * (1-psiI) * colSums(P * c(eta))
#' citaModI   <- round_preserve_sum(citaModI * N)
#' DBcitaModI <- createBase(citaModI)
#' 
#' # Creating auxiliary information
#' DBcitaModI[,AuxVar := rnorm(nrow(DBcitaModI), mean = 45, sd = 10)]
#' # Selects a sample with unequal probabilities
#' res <- S.piPS(n = 1200, as.data.frame(DBcitaModI)[,"AuxVar"])
#' sam <- res[,1]
#' pik <- res[,2]
#' DBcitaModISam <- copy(DBcitaModI[sam,])
#' DBcitaModISam[,Pik := pik]
#' 
#' # Gross flows estimation
#' estima <- estGF(sampleBase = DBcitaModISam, niter = 500, model = "II", colWeights = "Pik")
#' # gross flows variance estimation
#' varEstima <- reSamGF(sampleBase = DBcitaModISam, type = "Bootstrap", nRepBoot = 100,
#' 						model = "II", niter = 101,  colWeights = "Pik")
#' varEstima
#' @importFrom dplyr sample_n %>% select_ group_by_ summarize_all funs
#' @importFrom stats var
#' @references 
#' Efron, B. (1979), `Computers and the theory of statistics: Thinking the unthinkable', \emph{SIAM review} \bold{21}(4), pp. 460-480. \cr
#' Quenouille, M. H. (1949), `Problems in plane sampling', \emph{The Annals of Mathematical Statistics} pp. 355-375. \cr
#' Tukey, J. W. (1958), `Bias and confidence in not-quite large samples', \emph{Annals of Mathematical Statistics} \bold{29}, pp. 614.
#' @export reSamGF
reSamGF <- function(sampleBase = NULL, nRepBoot = 500, model = "I",
					niter = 100, type = "Bootstrap", colWeights = NULL,
					nonrft = FALSE){
	
	extractGF <- function(y) y[["Est.GF"]]
	extractRhoRR <- function(y) y[["Params.Model"]][["Est.rhoRR"]]
	varJack   <- function(y){
					temp <- (length(y) / (length(y)-1)) * sum((y - mean(y))^2)
				  	return(temp)
	}  


	if(!type %in% c("Bootstrap", "Jackknife")) stop("Your resampling method: '", type, "' is not supported")
	if(is.null(colWeights)) stop("Please specify a column for sampling weights")
	
	sampleBase 	<- as.data.frame(sampleBase)
	resulConsol <- NULL
	Nbase		<- nrow(sampleBase)
	
	if(type == "Bootstrap"){

		estBoot <- function(x){
			sampleBaseB  <- sample_n(sampleBase, size = Nbase, replace = TRUE)
 			r <- estGF(sampleBase = sampleBaseB, niter = niter, model = model, colWeights = colWeights, nonrft = nonrft)
					return(r)
		}
		it <- 1:nRepBoot
		resulConsol <- lapply(it, estBoot)
		
	}else if(type == "Jackknife"){

		estJack <- function(x){
			sampleBaseJ <- sampleBase[-x,]
			sampleBaseJ[,colWeights]  <- sampleBaseJ[,colWeights] + (sampleBase[x, colWeights] / (Nbase - 1))
			r <- estGF(sampleBase = sampleBaseJ, niter = niter, model = model, colWeights = colWeights, nonrft = nonrft)
			return(r)
		}
		it <- 1:Nbase
		resulConsol <- lapply(it, estJack)

	}
	
	acomSal <- function(x){
		x <- as.data.frame(x)
		x[,"id"] <- row.names(x)
		row.names(x) <- NULL
		return(x)
	}

	extractOP <- function(y,z){
		x <- as.data.frame(y[["Params.Model"]][[z]])
		x[,"id"] <- row.names(x)
		row.names(x) <- NULL
		names(x)[1] <- z
		return(x)
	}

	extractPij <- function(y){
		x <- as.data.frame(y[["Params.Model"]][["Est.pij"]])
		x[,"id"] <- row.names(x)
		row.names(x) <- NULL
		return(x)
	}
 
	varParam <- function(x, var1 = "Iteration", var2 = "id"){
		x[,var1] <- NULL
		y	<- x %>% 
			   	group_by_(var2) %>%
		  	   	summarize_all(ifelse(type == "Jackknife", funs(varJack), funs(var))) %>%
		  	   	as.data.frame()
		return(y) 					  
	}

	EstEta<- rbindlist(lapply(resulConsol,extractOP, "Est.eta"), idcol = "Iteration")
	EstPij<- rbindlist(lapply(resulConsol,extractPij), idcol = "Iteration")
	EstrhoRR <- rbindlist(lapply(lapply(resulConsol,extractRhoRR), acomSal), idcol = "Iteration")
	EstGF <- rbindlist(lapply(lapply(resulConsol,extractGF), acomSal), idcol = "Iteration")
	vGF  <- varParam(EstGF)
	vEta <- varParam(EstEta)
	vPij <- varParam(EstPij)
	vrhoRR <- varParam(EstrhoRR)

	if(!nonrft){
		EstrhoMM <- rbindlist(lapply(resulConsol,extractOP, "Est.rhoMM"), idcol = "Iteration")
		EstPsi   <- rbindlist(lapply(resulConsol,extractOP, "Est.psi"), idcol = "Iteration")
		vRhoMM  <- varParam(EstrhoMM)
		vPsi    <- varParam(EstPsi)
		return(list(vGF = vGF, vEta = vEta, vPij = vPij,
					vrhoRR = vrhoRR, vRhoMM = vRhoMM, vPsi = vPsi))
	}else{
		return(list(vGF = vGF, vEta = vEta, vPij = vPij, vrhoRR = vrhoRR))
	}

}
