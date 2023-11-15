#' Random generation from the pushed beta distribution
#'
#' \code{rpushbeta} returns random variables from the distribution.
#'
#' This function generates random variables from the pushed beta distribution.  Further details on the 
#' distribution can be found in the following paper:
#'
#' O'Neill, B. (2023) The pushed beta distribution and contaminated binary sampling.
#'
#' @usage \code{rpushbeta(n, shape1, shape2, shape3, push, left = TRUE)}
#' @param n The number of random variables to generate
#' @param shape1 The shape1 parameter for the pushed beta distribution (positive numeric value)
#' @param shape2 The shape2 parameter for the pushed beta distribution (positive numeric value)
#' @param shape3 The push-shape parameter for the pushed beta distribution (non-negative numeric value)
#' @param push The push parameter for the pushed beta distribution (numeric value between zero and one)
#' @param left Logical direction value; \code{TRUE} uses the left-pushed beta distribution; \code{FALSE} uses the right-pushed beta distribution
#' @return If all inputs are correctly specified (i.e., parameters are in allowable range and arguments are numeric)
#' then the output will be a vector of random variables from the distribution

rpushbeta <- function(n, shape1, shape2, shape3, push, left = TRUE) {

  #Check that inputs are appropriate type
  if (!is.numeric(n))                       stop('Error: Argument n is not numeric')
  if (!is.numeric(shape1))                  stop('Error: Shape1 parameter is not numeric')
  if (!is.numeric(shape2))                  stop('Error: Shape2 parameter is not numeric')
  if (!is.numeric(shape3))                  stop('Error: Push-shape parameter is not numeric')
  if (!is.numeric(push))                    stop('Error: Push parameter is not numeric')
  if (!is.logical(left))                    stop('Error: Direction input is not a logical value')
  
  #Check that parameters and options are atomic
  if (length(n) != 1)                       stop('Error: Argument n should be a single number')
  if (length(shape1) != 1)                  stop('Error: Shape1 parameter should be a single number')
  if (length(shape2) != 1)                  stop('Error: Shape2 parameter should be a single number')
  if (length(shape3) != 1)                  stop('Error: Shape3 (push-shape) parameter should be a single number')
  if (length(push)   != 1)                  stop('Error: Push parameter should be a single number')
  if (length(left)   != 1)                  stop('Error: Direction input should be a single logical value')
  
  #Check input n
  nn <- as.integer(n)
  if (nn != n)                              stop('Error: Argument n should be an integer')
  if (nn < 0)                               stop('Error: Argument n should be a non-negative integer')
  if (nn == Inf)                            stop('Error: Argument n should be finite')
  
  #Check that parameters are in allowable range
  if (min(shape1) <= 0)                     stop('Error: Shape1 parameter should be positive')
  if (min(shape2) <= 0)                     stop('Error: Shape2 parameter should be positive')
  if (min(shape3) < 0)                      stop('Error: Shape3 (push-shape) parameter should be non-negative')
  if (min(push) < 0)                        stop('Error: Push parameter should be between zero and one')
  if (max(push) > 1)                        stop('Error: Push parameter should be between zero and one')
  
  ######################################################################################################
  
  #Deal with the special case where distribution reduces to the beta distribution
  if (left) {
    if (shape3 == 0) { return(rbeta(n, shape1 = shape1, shape2 = shape2)) }
    if (push == 0)   { return(rbeta(n, shape1 = shape1, shape2 = shape2)) }
    if (push == 1)   { return(rbeta(n, shape1 = shape1, shape2 = shape2 + shape3)) } }
  if (!left) {
    if (shape3 == 0) { return(1-rbeta(n, shape1 = shape2, shape2 = shape1)) }
    if (push == 0)   { return(1-rbeta(n, shape1 = shape2, shape2 = shape1)) }
    if (push == 1)   { return(1-rbeta(n, shape1 = shape2 + shape3, shape2 = shape1)) } }
  
  #Generate random variables (left-pushed beta distribution)
  if (left) {
    
    #Compute scaling constant (using integration)
    KERNEL <- function(xx) { (xx^(shape1-1))*((1-xx)^(shape2-1))*((1-xx*push)^shape3) }
    INT    <- integrate(f = KERNEL, lower = 0, upper = 1)
    if (INT$message != "OK")   stop('Error: There was an error conducting numerical integration of the density kernel')
    if (INT$value <= 0)        stop('Error: There was an error conducting numerical integration of the density kernel')
    LOGSCALE <- log(INT$value)
    
    #Compute quantiles
    OUT   <- rep(0, n)
    LOG.P <- log(runif(n))
    LOG.P.KERN <- LOG.P + LOGSCALE
    for (i in 1:n) {
      
      #Compute quantiles from nonlinear optimisation
      kk <- LOG.P.KERN[i]
      q0 <- qbeta(LOG.P[i], shape1 = shape1, shape2 = shape2 + push*shape3, lower.tail = TRUE, log.p = TRUE)
      x0 <- log(q0/(1-q0))
      OBJECTIVE <- function(xx) {
        qq  <- exp(xx)/(1+exp(xx))
        INT <- integrate(f = KERNEL, lower = 0, upper = qq)
        (log(INT$value) - kk)^2 }
      NLM <- nlm(f = OBJECTIVE, p = x0, gradtol = 10e-21)
      XX <- NLM$estimate
      OUT[i] <- exp(XX)/(1+exp(XX)) } }

  #Compute quantile vector (right-pushed beta distribution)
  if (!left) {
    
    #Compute scaling constant (using integration)
    KERNEL <- function(xx) { (xx^(shape1-1))*((1-xx)^(shape2-1))*((1-push+xx*push)^shape3) }
    INT    <- integrate(f = KERNEL, lower = 0, upper = 1)
    if (INT$message != "OK")   stop('Error: There was an error conducting numerical integration of the density kernel')
    if (INT$value <= 0)        stop('Error: There was an error conducting numerical integration of the density kernel')
    LOGSCALE <- log(INT$value)
    
    #Compute quantiles
    OUT   <- rep(0, n)
    LOG.P <- log(runif(n))
    LOG.P.KERN <- LOG.P + LOGSCALE
    for (i in 1:n) {
      
      #Compute quantiles from nonlinear optimisation
      kk <- LOG.P.KERN[i]
      q0 <- qbeta(LOG.P[i], shape1 = shape1, shape2 = shape2 + push*shape3, lower.tail = TRUE, log.p = TRUE)
      x0 <- log(q0/(1-q0))
      OBJECTIVE <- function(xx) {
        qq  <- exp(xx)/(1+exp(xx))
        INT <- integrate(f = KERNEL, lower = 0, upper = qq)
        (log(INT$value) - kk)^2 }
      NLM <- nlm(f = OBJECTIVE, p = x0, gradtol = 10e-21)
      XX <- NLM$estimate
      OUT[i] <- exp(XX)/(1+exp(XX)) } }

  #Return output
  OUT }
