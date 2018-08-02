#' Gross Flows estimation
#' @description Gross Flows under complex electoral surveys.
#' @param sampleBase An object of class "data.frame" containing the information of electoral candidates. The data must contain the samplings weights.
#' @param niter The number of iterations for the \eqn{\eta_{i}} and \eqn{p_{ij}} model parameters within the model.
#' @param colWeights The column name containing the sampling weights to be used in the fitting process.
#' @param model A character indicating the model to be used in estimating estimated gross flows. The models available are: "I","II","III","IV" (see also "Details").
#' @param nonrft A logical value indicating a non response for first time. 
#' @return \code{estGF} returns a list containing:
#' 			\enumerate{
#' 			 \item \bold{Est.CIV:} a data.frame containing the gross flows estimation.
#' 			 \item \bold{Params.Model:} a list that contains the \eqn{\hat{\eta}_{i}}, \eqn{\hat{p}_{ij}}, \eqn{\hat{\psi}(i,j)}, \eqn{\hat{\rho}_{RR}(i,j)}, \eqn{\hat{\rho}_{MM}(i,j)} parameters for the estimated model.
#' 			 \item \bold{Sam.Est:} a list containing the sampling estimators \eqn{\hat{N}_{ij}}, \eqn{\hat{R}_{i}}, \eqn{\hat{C}_{j}}, \eqn{\hat{M}}, \eqn{\hat{N}}.
#' 			}
#' @details
#' The population size \eqn{N} must satisfy the condition:
#' \deqn{ N = \sum_{j}\sum_{i} N_{ij} + \sum_{j} C_{j} + \sum_{i} R_{i} + M}
#' where, \eqn{N_{ij}} is the amount of people interviewed who have classification \eqn{i} at first time and classification \eqn{j} at second time, \eqn{R_{i}} is the amount of people who did not respond at second time, but did at first time, \eqn{C_{j}} is the amount of people who did not respond at first time, but they did at second time and \eqn{M} is the number of people who did not respond at any time or could not be reached.
#' Let \eqn{\eta_{i}} the initial probability that a person has classification \eqn{i} in the first time, and let \eqn{p_{ij}} the vote transition probability for the cell \eqn{(i,j)}, where \eqn{\sum_{i} \eta_{i} = 1} and \eqn{\sum_{j} p_{ij} = 1}.
#' Thus, four possibles models for the gross flows are given by:
#' \enumerate{
#' 	\item \bold{Model I:} This model assumes that a person's initial probability of being classified as \eqn{i} at first time is the same for everyone, that is, \eqn{\psi(i,j) = \psi}. Besides, transition probabilities between respond and non response not depend of the classification \eqn{(i,j)}, that is \eqn{\rho_{MM}(i,j) = \rho_{MM}} and \eqn{\rho_{RR}(i,j) = \rho_{RR}}.
#'  \item \bold{Model II:} Unlike 'Model I', this model assumes that person initial probability that person has classification \eqn{(i,j)}, only depends of his classification at first time, that is \eqn{\psi(i,j) = \psi(i)}.
#' 	\item \bold{Model III:} Unlike 'Model I', this model assumes that transition probabilities between response and non response only depends of probability classification at first time, that is \eqn{\rho_{MM}(i,j) = \rho_{MM}(i)} and \eqn{\rho_{RR}(i,j) = \rho_{RR}(i)}.
#' 	\item \bold{Model IV:} Unlike 'Model I', this model assumes that transition probabilities between response and non response only depends of probability classification at second time, that is \eqn{\rho_{MM}(i,j) = \rho_{MM}(j)} and \eqn{\rho_{RR}(i,j) = \rho_{RR}(j)}.
#' }
#' @examples
#' library(TeachingSampling)
#' library(data.table)
#' # Colombia's electoral candidates in 2014
#' candidates_t0 <- c("Clara","Enrique","Santos","Martha","Zuluaga","WhiteVote", "NoVote")
#' candidates_t1 <- c("Santos","Zuluaga","WhiteVote", "NoVote")
#' 
#' N <- 100000
#' nCanT0 <- length(candidates_t0)
#' nCanT1 <- length(candidates_t1)
#' # Initial probabilities
#' eta <- matrix(c(0.10, 0.10, 0.20, 0.17, 0.28, 0.1, 0.05),
#' 				byrow = TRUE, nrow = nCanT0)
#' # Transition probabilities
#' P <- matrix(c(0.10, 0.60, 0.15, 0.15,
#' 				 0.30, 0.10, 0.25,0.35,
#' 				 0.34, 0.25, 0.16, 0.25,
#' 				 0.25,0.05, 0.35,0.35,
#' 				 0.10, 0.25, 0.45,0.20,
#' 				 0.12, 0.36, 0.22, 0.30,
#' 				 0.10,0.15, 0.30,0.45),
#' 		byrow = TRUE, nrow = nCanT0)
#' citaMod <- matrix(, ncol = nCanT1, nrow = nCanT0)
#' row.names(citaMod) <- candidates_t0
#' colnames(citaMod) <- candidates_t1
#' 
#' for(ii in 1:nCanT0){
#' 		citaMod[ii,] <- c(rmultinom(1, size = N * eta[ii,], prob = P[ii,]))
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
#' citaModI[1:nCanT0, 1:nCanT1] <- P * c(eta) * rhoRRI * psiI
#' citaModI[(nCanT0 + 1), (nCanT1 + 1)] <- rhoMMI * (1-psiI)
#' citaModI[1:nCanT0, (nCanT1 + 1)] <- (1-rhoRRI) * psiI * rowSums(P * c(eta))
#' citaModI[(nCanT0 + 1), 1:nCanT1 ] <- (1-rhoMMI) * (1-psiI) * colSums(P * c(eta))
#' citaModI <- round_preserve_sum(citaModI * N)
#' DBcitaModI <- createBase(citaModI)
#' 
#' # Creating auxiliary information
#' DBcitaModI[,AuxVar := rnorm(nrow(DBcitaModI), mean = 45, sd = 10)]
#' 
#' # Selects a sample with unequal probabilities
#' res <- S.piPS(n = 3200, as.data.frame(DBcitaModI)[,"AuxVar"])
#' sam <- res[,1]
#' pik <- res[,2]
#' DBcitaModISam <- copy(DBcitaModI[sam,])
#' DBcitaModISam[,Pik := pik]
#' 
#' # Gross Flows estimation
#' estima <- estGF(sampleBase = DBcitaModISam, niter = 500, model = "I", colWeights = "Pik")
#' estima
#' @references 
#' Stasny, E. (1987), `Some markov-chain models for nonresponse in estimating gross', \emph{Journal of Oficial Statistics} \bold{3}, pp. 359-373. \cr
#' Sarndal, C.-E., Swensson, B. \& Wretman, J. (1992), \emph{Model Assisted Survey Sampling}, Springer-Verlag, New York, USA. \cr
#' Gutierrez, A., Trujillo, L. \& Silva, N. (2014), `The estimation of gross ows in complex surveys with random nonresponse', \emph{Survey Methodology} \bold{40}(2), pp. 285-321.
#' @export estGF
estGF  <- function(sampleBase = NULL, niter = 100, model = NULL, colWeights = NULL, nonrft = FALSE){
	
	if(is.null(sampleBase)) stop("Data base must be specified")
	if(niter == 100) warning("By defect ", niter, " iterations will be run for: Eta and Pij", call. = FALSE)
	if(!model %in% c('I', 'II', 'III', 'IV')) stop("The model is not supported by method")
	if(is.null(colWeights)) stop("Please specify a column for sampling weights")
	if(!class(sampleBase)[1] %in% c("data.frame", "data.table")) stop("str dataBase must be data.frame or data.table")
	if(nonrft) warning('nonresponse at first time is ommited')

	sampleBase 	  <- as.data.frame(sampleBase)
	candidates_t0 <- grep("^(t0).*",names(sampleBase), value = TRUE)
	candidates_t0 <- candidates_t0[!grepl("t0_Non_Resp", candidates_t0)]
	candidates_t1 <- grep("^(t1).*",names(sampleBase), value = TRUE)
	candidates_t1 <- candidates_t1[!grepl("t1_Non_Resp", candidates_t1)]
	NCandidatosT0 <- length(candidates_t0)
	NCandidatosT1 <- length(candidates_t1)

	if(!"t0_Non_Resp" %in% names(sampleBase)) sampleBase[,"t0_Non_Resp"] <- 0
	if(!"t1_Non_Resp" %in% names(sampleBase)) sampleBase[,"t1_Non_Resp"] <- 0

	z1 <- (1 - sampleBase[,"t0_Non_Resp"])
	z2 <- (1 - sampleBase[,"t1_Non_Resp"])
	y1 <- as.matrix(sampleBase[,candidates_t0])
	y2 <- as.matrix(sampleBase[,candidates_t1])
	Wk <- 1/sampleBase[,colWeights]
		
	Est.Nij <- t(Wk*y1)%*%y2 
	Est.Ri  <- t(Wk*y1)%*%(1-z2) 
	Est.Cj  <- t(Wk*y2)%*%(1-z1) 
	Est.M   <- c(t(Wk*(1-z1))%*%(1-z2) )
	Est.N   <- sum(Est.Nij)+sum(Est.Ri)+sum(Est.Cj)+Est.M 

	if(model == "I" ){
	
		Est.rhoRR <- sum(Est.Nij)/(sum(Est.Nij)+sum(Est.Ri))
		Est.pij   <- as.matrix(prop.table(Est.Nij, 1))
		Est.rhoMM <- NA
		Est.psi   <- 1		
		Est.eta   <- (rowSums(Est.Nij) + Est.Ri) / Est.N
		
		if(!nonrft){
			Est.psi   <- (sum(Est.Nij)+sum(Est.Ri))/(sum(Est.Nij)+sum(Est.Ri)+sum(Est.Cj)+Est.M)
			Est.rhoMM <- Est.M/(sum(Est.Cj)+Est.M)
			Est.eta <- rowSums(Est.Nij) / sum(Est.Nij)
		
			for(kk in 1:niter){
		
				M.eta     <- NULL
				for(ii in 1:NCandidatosT0){
					etai <- sum(t(Est.Cj) * Est.pij[ii,] / (Est.eta %*% Est.pij))
					M.eta <- c(M.eta, etai)				
				}
				
				Est.eta.aux <- (rowSums(Est.Nij) + Est.Ri + 
								matrix(M.eta * Est.eta)) / (sum(Est.Nij) + 
								sum(Est.Ri)+sum(Est.Cj) )
		
				parte1  <- t(kronecker(t(matrix(Est.eta)), matrix(t(Est.Cj)/(Est.eta%*%Est.pij))))
				parte2  <- matrix(rep(rowSums(Est.Nij),NCandidatosT1), nrow = NCandidatosT0)
				parte3  <- matrix(rep(M.eta*Est.eta,NCandidatosT1), nrow = NCandidatosT0)
				Est.pij <- (Est.Nij+Est.pij*parte1)/(parte2+parte3)
				Est.eta <- t(Est.eta.aux)

			}
			Est.eta <- t(Est.eta)
		}
			
		Est.GF <- as.vector(sum(Est.Nij))*as.matrix(Est.pij)*as.vector(Est.eta)

		}else if(model == "II"){

			Est.rhoRR <- sum(Est.Nij)/(sum(Est.Nij)+sum(Est.Ri))
			Est.pij   <- as.matrix(prop.table(Est.Nij, 1))
			Est.rhoMM <- NA
			Est.psi   <- 1	
			Est.eta   <- (rowSums(Est.Nij) + Est.Ri) / Est.N
			
			if(!nonrft){
				Est.rhoMM <- Est.M/(sum(Est.Cj)+Est.M)
				Est.psi   <- rep(1/NCandidatosT0, NCandidatosT0)
				Est.eta   <- rowSums(Est.Nij) / sum(Est.Nij)
						
				for(kk in 1:niter){

					suma.psi <- NULL
					suma.eta <- NULL
					for(jj in 1:NCandidatosT1){
						suma.psiAux <- Est.Cj[jj]*(1-Est.psi)*Est.eta*Est.pij[,jj]/sum((1-Est.psi)*Est.eta*Est.pij[,jj])
						suma.psi 	<- cbind(suma.psi, suma.psiAux)
						
					}

					suma.psi <- rowSums(suma.psi)
					Est.psi  <- as.vector((rowSums(Est.Nij) + Est.Ri) / (rowSums(Est.Nij) + 
											Est.Ri + suma.psi + Est.M*((1-Est.psi)*Est.eta /
											 sum((1-Est.psi)*Est.eta))))

					for(jj in 1:NCandidatosT1){
					suma.etaAux	<- Est.Cj[jj]*(1-Est.psi)*Est.eta*Est.pij[,jj]/sum((1-Est.psi)*Est.eta*Est.pij[,jj])
					suma.eta 	<- cbind(suma.eta, suma.etaAux)
					}

					suma.eta <- rowSums(suma.eta)
					Est.eta <- as.vector((rowSums(Est.Nij) + c(Est.Ri) + suma.eta + 
										    Est.M*((1-Est.psi)*Est.eta / sum((1-Est.psi)*Est.eta))) /
											  (sum(Est.Nij)+sum(Est.Ri)+sum(Est.Cj)+Est.M))

					M.pij <- NULL
					for(jj in 1:NCandidatosT1){
						M.pijAux <- sum((1-Est.psi)*Est.eta*Est.pij[,jj])
						M.pij 	 <- c(M.pij, M.pijAux)
					}

					matriz1.pij <- kronecker(t(matrix(Est.Cj)), (1-Est.psi)*Est.eta)*Est.pij
					matriz2.pij <- kronecker(t(matrix(M.pij)),
											 matrix(rep(1,NCandidatosT0)))
					matriz3.pij <- kronecker(t(matrix(rep(1,NCandidatosT1))), rowSums(Est.Nij+(matriz1.pij/matriz2.pij)))
					Est.pij <- (Est.Nij+(matriz1.pij/matriz2.pij))/matriz3.pij
				}
			}

			Est.GF <- as.vector(sum(Est.Nij))*as.matrix(Est.pij)*as.vector(Est.eta)
			names(Est.eta) <- candidates_t0
			Est.eta <- t(t(Est.eta))
			if(!nonrft){
				names(Est.psi) <- candidates_t0
				Est.psi <- t(t(Est.psi))
			}
			
		}else if(model == "III"){

				Est.rhoRR <- rowSums(Est.Nij)/(rowSums(Est.Nij)+Est.Ri)
				Est.pij   <- as.matrix(prop.table(Est.Nij, 1))
				Est.eta   <- (rowSums(Est.Nij) + Est.Ri) / Est.N
				Est.rhoMM <- NA
				Est.psi   <- 1 
				
				if(!nonrft){
					Est.rhoMM <- rep(1/NCandidatosT0, NCandidatosT0)
					Est.psi   <- (sum(Est.Nij)+sum(Est.Ri))/(sum(Est.Nij)+sum(Est.Ri)+sum(Est.Cj)+Est.M)
					Est.eta   <- rowSums(Est.Nij) / sum(Est.Nij)
								
					for(kk in 1:niter){

						suma.rhoMM <- NULL
						for(ii in 1:NCandidatosT1){
							suma.rhoMMAux <- Est.Cj[ii]*(1-Est.rhoMM)*Est.eta*Est.pij[,ii]/sum((1-Est.rhoMM)*Est.eta*Est.pij[,ii])
							suma.rhoMM   <- cbind(suma.rhoMM, suma.rhoMMAux)
						}

						P1.rhoMM <- rowSums(suma.rhoMM)
						P2.rhoMM <- Est.M*Est.rhoMM*Est.eta/sum(Est.rhoMM*Est.eta)
						Est.rhoMM <- P2.rhoMM/(P1.rhoMM+P2.rhoMM)
						if(sum(is.nan(Est.rhoMM)) == NCandidatosT0) Est.rhoMM[is.nan(Est.rhoMM)] <- 0 # Pilas con esto

						suma.eta <- NULL
						for(jj in 1:NCandidatosT1){
							suma.etaAux	<- Est.Cj[jj]*(1-Est.rhoMM)*Est.eta*Est.pij[,jj]/sum((1-Est.rhoMM)*Est.eta*Est.pij[,jj])
							suma.eta 	<- cbind(suma.eta, suma.etaAux)
						}

						P1.eta <- rowSums(suma.rhoMM)
						P2.eta <- Est.M*Est.rhoMM*Est.eta/sum(Est.rhoMM*Est.eta)
						Est.eta <- (rowSums(Est.Nij)+c(Est.Ri)+P1.eta+P2.eta)/(sum(Est.Nij)+sum(Est.Ri)+sum(Est.Cj)+Est.M)

						M.pij <- NULL
						for(jj in 1:NCandidatosT1){
							M.pijAux <- sum((1-Est.psi)*Est.eta*Est.pij[,jj])
							M.pij 	 <- c(M.pij, M.pijAux)
						}

						matriz1.pij <- kronecker(t(matrix(Est.Cj)), (1-Est.rhoRR)*Est.eta)*Est.pij
						matriz2.pij <- kronecker(t(matrix(M.pij)),matrix(rep(1,NCandidatosT0)))
						matriz3.pij <- kronecker(t(rep(1,NCandidatosT1)), rowSums(Est.Nij+(matriz1.pij/matriz2.pij)))
						Est.pij <- (Est.Nij+(matriz1.pij/matriz2.pij))/matriz3.pij
					}
				}
				Est.rhoMM <- t(t(Est.rhoMM))
				Est.eta <- t(t(Est.eta))
				Est.GF <- as.vector(sum(Est.Nij))*as.matrix(Est.pij)*as.vector(Est.eta)		

		}else if(model == "IV"){

			if(nonrft){
				Est.rhoMM <- NA
				Est.psi   <- 1
				Est.eta   <- (rowSums(Est.Nij) + Est.Ri) / Est.N
			}else{
				Est.psi <- (sum(Est.Nij)+sum(Est.Ri))/(sum(Est.Nij)+sum(Est.Ri)+sum(Est.Cj)+Est.M)
				Est.rhoMM <- rep(1/NCandidatosT0, NCandidatosT1)
				Est.eta   <- rowSums(Est.Nij) / sum(Est.Nij)
			}

				Est.rhoRR <- rep(1/NCandidatosT0, NCandidatosT1)
				Est.pij   <- as.matrix(prop.table(Est.Nij, 1))

				for(kk in 1:niter){			

					suma.rhoRR <- NULL
					for(ii in 1:NCandidatosT0){
						suma.rhoRRAux <- Est.Ri[ii]*(1-Est.rhoRR)*Est.pij[ii,]/sum((1-Est.rhoRR)*Est.pij[ii,])
						suma.rhoRR   <- rbind(suma.rhoRR, suma.rhoRRAux)
					}

					P1.rhoRR <- colSums(suma.rhoRR)
					Est.rhoRR <- colSums(Est.Nij)/(colSums(Est.Nij)+P1.rhoRR)

					if(!nonrft){
						Est.rhoMM <- 1-(Est.Cj*sum(kronecker(t(Est.rhoMM), Est.eta)*Est.pij)/(Est.M*colSums(Est.eta*Est.pij)))
						Est.rhoMM <- c(Est.rhoMM)

						suma.eta <- NULL
						for(jj in 1:NCandidatosT1){
							suma.etaAux	<- Est.Cj[jj]*Est.eta*Est.pij[,jj]/sum(Est.eta*Est.pij[,jj])
							suma.eta 	<- cbind(suma.eta, suma.etaAux)
						}

						P1.eta <- rowSums(suma.eta)
						P2.eta <- Est.M*rowSums(t(Est.rhoMM*t(Est.eta*Est.pij)))/sum(Est.rhoMM*t(Est.eta*Est.pij))
						Est.eta <- (rowSums(Est.Nij)+c(Est.Ri)+P1.eta+P2.eta)/(sum(Est.Nij)+sum(Est.Ri)+sum(Est.Cj)+Est.M)
					}

					matriz1.R.pij <- kronecker(t((1-Est.rhoRR)), matrix(Est.Ri))*Est.pij
					M.R.pij <- NULL
					for(ii in 1:NCandidatosT0){
						M.R.pijAux <- sum((1-Est.rhoRR)*Est.pij[ii,])
						M.R.pij    <- c(M.R.pij, M.R.pijAux)
					}

					matriz2.R.pij  <- kronecker(t(matrix(rep(1,NCandidatosT1))), matrix(M.R.pij) )
					P1 <- matriz1.R.pij/matriz2.R.pij

					if(!nonrft){		
						matriz1.C.pij  <- kronecker(t(matrix(Est.Cj)), Est.eta)*Est.pij
						M.C.pij <- NULL
						for(jj in 1:NCandidatosT1){
							M.C.pijAux <- sum(Est.eta*Est.pij[,jj])
							M.C.pij    <- c(M.C.pij, M.C.pijAux)
						}
						matriz2.C.pij  <- kronecker(t(matrix(M.C.pij)), matrix(rep(1,NCandidatosT0)))
						P2 <-matriz1.C.pij/matriz2.C.pij
						P3 <- c(Est.M)*t(Est.rhoMM*t(Est.eta*Est.pij))/sum(t(Est.rhoMM*t(Est.eta*Est.pij)))
						P4 <- rowSums(Est.Nij+P1+P2+P3)
						P5 <- kronecker(t(matrix(rep(1,NCandidatosT1))), P4)
						Est.pij <- (Est.Nij+P1+P2+P3)/P5
					}else{
					 	Est.pij <- (Est.Nij+P1) / c(rowSums(Est.Nij) + Est.Ri)
					}
				}
				Est.rhoMM <- t(Est.rhoMM)
				Est.rhoRR <- t(Est.rhoRR)
				Est.eta   <- t(t(Est.eta))
				Est.GF <- as.vector(sum(Est.Nij))*as.matrix(Est.pij)*as.vector(Est.eta)
				if(!nonrft){
					names(Est.rhoMM) <- candidates_t1
					Est.rhoMM <- t(Est.rhoMM)
				}	
		} 

		if(!nonrft){
			return(list(Est.GF = Est.GF,
						Params.Model = list(Est.eta = Est.eta, Est.pij = Est.pij,
									        Est.psi = Est.psi, Est.rhoRR = Est.rhoRR,
									        Est.rhoMM = Est.rhoMM),
						Sam.Est  =  list(Est.Nij = Est.Nij, Est.Ri = Est.Ri,
										 Est.Cj = Est.Cj, Est.M = Est.M,
										 Est.N = Est.N)))
		}else{
			return(list(Est.GF = Est.GF,
						Params.Model = list(Est.eta = Est.eta, Est.pij = Est.pij,
											Est.rhoRR = Est.rhoRR),
						Sam.Est  =  list(Est.Nij = Est.Nij, Est.Ri = Est.Ri,
										 Est.N = Est.N)))
		}

}