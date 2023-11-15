#' Density function of the pushed beta distribution
#'
#' \code{dpushbeta} returns the density or log-density values for the arguments.
#'
#' This function computes densities or log-densities from the density function of the pushed beta
#' distribution.  Further details on the distribution can be found in the following paper:
#'
#' O'Neill, B. (2023) The pushed beta distribution and contaminated binary sampling.
#'
#' @usage \code{dpushbeta(x, shape1, shape2, shape3, push, left = TRUE, log = FALSE)}
#' @param x A vector of numeric values to be used as arguments for the density function
#' @param shape1 The shape1 parameter for the pushed beta distribution (positive numeric value)
#' @param shape2 The shape2 parameter for the pushed beta distribution (positive numeric value)
#' @param shape3 The push-shape parameter for the pushed beta distribution (non-negative numeric value)
#' @param push The push parameter for the pushed beta distribution (numeric value between zero and one)
#' @param left Logical direction value; \code{TRUE} uses the left-pushed beta distribution; \code{FALSE} uses the right-pushed beta distribution
#' @param log A logical value specifying whether results should be returned as log-densities
#' @return If all inputs are correctly specified (i.e., parameters are in allowable range and arguments are numeric)
#' then the output will be a vector of densities/log-densities corresponding to the vector argument x

dpushbeta <- function(x, shape1, shape2, shape3, push, left = TRUE, log = FALSE) {

  #Check that inputs are appropriate type
  if (!is.numeric(x))                       stop('Error: Argument x is not numeric')
  if (!is.numeric(shape1))                  stop('Error: Shape1 parameter is not numeric')
  if (!is.numeric(shape2))                  stop('Error: Shape2 parameter is not numeric')
  if (!is.numeric(shape3))                  stop('Error: Push-shape parameter is not numeric')
  if (!is.numeric(push))                    stop('Error: Push parameter is not numeric')
  if (!is.logical(left))                    stop('Error: Direction input (left) is not a logical value')
  if (!is.logical(log))                     stop('Error: log option is not a logical value')
  
  #Check that parameters and options are atomic
  if (length(shape1) != 1)                  stop('Error: Shape1 parameter should be a single number')
  if (length(shape2) != 1)                  stop('Error: Shape2 parameter should be a single number')
  if (length(shape3) != 1)                  stop('Error: Shape3 (push-shape) parameter should be a single number')
  if (length(push)   != 1)                  stop('Error: Push parameter should be a single number')
  if (length(left)   != 1)                  stop('Error: Direction input (left) should be a single logical value')
  if (length(log)    != 1)                  stop('Error: log option should be a single logical value')

  #Check that parameters are in allowable range
  if (min(shape1) <= 0)                     stop('Error: Shape1 parameter should be positive')
  if (min(shape2) <= 0)                     stop('Error: Shape2 parameter should be positive')
  if (min(shape3) < 0)                      stop('Error: Shape3 (push-shape) parameter should be non-negative')
  if (min(push) < 0)                        stop('Error: Push parameter should be between zero and one')
  if (max(push) > 1)                        stop('Error: Push parameter should be between zero and one')
  
  ######################################################################################################
  
  #Deal with the special case where distribution reduces to the beta distribution
  if (left) {
    if (shape3 == 0) { return(dbeta(x, shape1 = shape1, shape2 = shape2, log = log)) }
    if (push == 0)   { return(dbeta(x, shape1 = shape1, shape2 = shape2, log = log)) }
    if (push == 1)   { return(dbeta(x, shape1 = shape1, shape2 = shape2 + shape3, log = log)) } }
  if (!left) {
    if (shape3 == 0) { return(dbeta(1-x, shape1 = shape2, shape2 = shape1, log = log)) }
    if (push == 0)   { return(dbeta(1-x, shape1 = shape2, shape2 = shape1, log = log)) }
    if (push == 1)   { return(dbeta(1-x, shape1 = shape2 + shape3, shape2 = shape1, log = log)) } }
  
  #Compute log-density vector (left-pushed beta distribution)
  if (left) {
    
    #Compute log-kernel values
    LOGKERNEL <- rep(-Inf, length(x))
    for (i in 1:length(x)) {
      xx <- x[i]
      if ((xx == 0)&(shape1 == 1)) { LOGKERNEL[i] <- 0 }
      if ((xx == 1)&(shape2 == 1)) { LOGKERNEL[i] <- shape3*log(1-push) }
      if ((xx > 0)&(xx < 1)) { 
        LOGKERNEL[i] <- (shape1-1)*log(xx) + (shape2-1)*log(1-xx) + shape3*log(1-xx*push) } }
    
    #Compute scaling constant (using integration)
    KERNEL <- function(xx) { (xx^(shape1-1))*((1-xx)^(shape2-1))*((1-xx*push)^shape3) }
    INT    <- integrate(f = KERNEL, lower = 0, upper = 1)
    if (INT$message != "OK")   stop('Error: There was an error conducting numerical integration of the density kernel')
    if (INT$value <= 0)        stop('Error: There was an error conducting numerical integration of the density kernel')
    LOGSCALE <- log(INT$value)
    
    #Generate log-density vector
    LOGDENSITY <- LOGKERNEL - LOGSCALE }
  
  #Compute log-density vector (right-pushed beta distribution)
  if (!left) {
    
    #Compute log-kernel values
    LOGKERNEL <- rep(-Inf, length(x))
    for (i in 1:length(x)) {
      xx <- x[i]
      if ((xx == 0)&(shape1 == 1)) { LOGKERNEL[i] <- shape3*log(1-push) }
      if ((xx == 1)&(shape2 == 1)) { LOGKERNEL[i] <- 0 }
      if ((xx > 0)&(xx < 1)) { 
        LOGKERNEL[i] <- (shape1-1)*log(xx) + (shape2-1)*log(1-xx) + shape3*log(1-push+xx*push) } }
    
    #Compute scaling constant (using integration)
    KERNEL <- function(xx) { (xx^(shape1-1))*((1-xx)^(shape2-1))*((1-push+xx*push)^shape3) }
    INT    <- integrate(f = KERNEL, lower = 0, upper = 1)
    if (INT$message != "OK")   stop('Error: There was an error conducting numerical integration of the density kernel')
    if (INT$value <= 0)        stop('Error: There was an error conducting numerical integration of the density kernel')
    LOGSCALE <- log(INT$value)
    
    #Generate log-density vector
    LOGDENSITY <- LOGKERNEL - LOGSCALE }

  #Return output
  if (log) { LOGDENSITY } else { exp(LOGDENSITY) } }

