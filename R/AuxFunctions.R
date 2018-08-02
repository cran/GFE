#' Create a database for gross flows.
#' @description Create a database based on \eqn{\xi} model.
#' @param x A matrix that contains information of the observable process.
#' @return \code{createBase} returns \code{data.table,data.frame} that contains the data                                                                                  base based on \eqn{\xi} model.
#' @examples
#' candidates_t0 <- c("Candidate1","Candidate2","Candidate3","Candidate4",
#' 	"Candidate5","WhiteVote", "NoVote")
#' candidates_t1 <- c("Candidate3","Candidate5","WhiteVote", "NoVote")
#' 
#' N 	   <- 100000
#' nCanT0  <- length(candidates_t0)
#' nCanT1  <- length(candidates_t1)
#' 
#' eta <- matrix(c(0.10, 0.10, 0.20, 0.17, 0.28, 0.1, 0.05),
#' 			   byrow = TRUE, nrow = nCanT0)
#' P <- matrix(c(0.10, 0.60, 0.15, 0.15,
#' 			  0.30, 0.10, 0.25,	0.35,
#' 			  0.34, 0.25, 0.16, 0.25,
#' 			  0.25,	0.05, 0.35,	0.35,
#' 			  0.10, 0.25, 0.45,	0.20,
#' 			  0.12, 0.36, 0.22, 0.30,
#' 			  0.10,	0.15, 0.30,	0.45),
#' 	 byrow = TRUE, nrow = nCanT0)
#' 
#' citaModel <- matrix(, ncol = nCanT1, nrow = nCanT0)
#' row.names(citaModel) <- candidates_t0
#' colnames(citaModel) <- candidates_t1
#' 
#' for(ii in 1:nCanT0){ 
#'  citaModel[ii,] <- c(rmultinom(1, size = N * eta[ii,], prob = P[ii,]))
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
#' DBmodCitaI <- createBase(citaModI)
#' DBmodCitaI
#' @import data.table
#' @importFrom TeachingSampling Domains
#' @export createBase
createBase <- function(x){
	candT0 <- rownames(x)
	candT1 <- rownames(t(x))
	t0  <- rep(rownames(x), apply(x, 1, sum))
	t0D <- data.table(Domains(t0))
	setcolorder(t0D, candT0)
		funRep <- function(x){
			aux <- rep(names(x), x)
			doAux <- data.table(Domains(aux))
		}
	t1D  <- rbindlist(apply(x, 1, funRep), fill = TRUE)
	setcolorder(t1D, candT1)
	colnames(t0D) <- paste("t0", candT0, sep ="_")
	colnames(t1D) <- paste("t1", candT1, sep ="_")
	base <- cbind(t0D, t1D)
	base[is.na(base),] <- 0
return(base)
}
#' Round preserve sum.
#' 
#' @description Rounds a vector of numbers while preserving the sum of them.
#' @param x A numeric vector.
#' @param digits The number of digits to take in account in the rounding process.
#' @return \code{round_preserve_sum} returns \code{y} with round vector.
#' @examples
#' sum(c(0.333, 0.333, 0.334))
#' round(c(0.333, 0.333, 0.334), 2)
#' sum(round(c(0.333, 0.333, 0.334), 2))
#' round_preserve_sum(c(0.333, 0.333, 0.334), 2)
#' sum(round_preserve_sum(c(0.333, 0.333, 0.334), 2))
#' @importFrom utils tail
#' @source
#' \url{https://www.r-bloggers.com/round-values-while-preserve-their-rounded-sum-in-r/} and \url{http://stackoverflow.com/questions/32544646/round-vector-of-numerics-to-integer-while-preserving-their-sum}
#' @export round_preserve_sum
round_preserve_sum <- function(x, digits = 0) {
  up <- 10 ^ digits
  x <- x * up
  y <- floor(x)
  indices <- tail(order(x-y), round(sum(x)) - sum(y))
  y[indices] <- y[indices] + 1
  return(y / up)
}
