#' @title Method of moments estimation ( beta distribution )
#' @usage beta.mom(qs.in)
#' @param qs.in A vector contains the numbers that will be fitted with a beta distribution.
#' @details beta.mom() function can be used to estimate parameters in a Beta function using method of moments
#' @return alpha.hat,beta.hat: Returns the estimation of alpha and beta.
#' @author Ning Leng
#' @examples beta.mom(rbeta(10,1,1))

beta.mom <-
function(qs.in){
	xbar <- mean(qs.in)
	s2 <- var(qs.in)
	term <- (xbar*(1-xbar))/s2
	alpha.hat <- xbar*(term-1)
	beta.hat <- (1-xbar)*(term-1)
	return(c(alpha.hat,beta.hat))
}

