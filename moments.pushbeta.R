#' Moments of the pushed beta distribution
#'
#' \code{moments.pushbeta} returns some representative moments from the distribution.
#'
#' This function computes some representative moments from the pushed beta distribution.  Further details 
#' on the distribution can be found in the following paper:
#'
#' O'Neill, B. (2023) The pushed beta distribution and contaminated binary sampling.
#'
#' @usage \code{moments.pushbeta(x, shape1, shape2, shape3, push, left = TRUE, lower.tail = TRUE, log.p = FALSE)}
#' @param shape1 The shape1 parameter for the pushed beta distribution (positive numeric value)
#' @param shape2 The shape2 parameter for the pushed beta distribution (positive numeric value)
#' @param shape3 The push-shape parameter for the pushed beta distribution (non-negative numeric value)
#' @param push The push parameter for the pushed beta distribution (numeric value between zero and one)
#' @param left Logical direction value; \code{TRUE} uses the left-pushed beta distribution; \code{FALSE} uses the right-pushed beta distribution
#' @param include.sd Logical value; if \code{TRUE} the output includes the standard deviation
#' @return If all inputs are correctly specified (i.e., parameters are in allowable range) then the output will be
#' a data frame of moments

moments.pushbeta <- function(shape1, shape2, shape3, push, left = TRUE, include.sd = FALSE) {

  #Check that inputs are appropriate type
  if (!is.numeric(shape1))                  stop('Error: Shape1 parameter is not numeric')
  if (!is.numeric(shape2))                  stop('Error: Shape2 parameter is not numeric')
  if (!is.numeric(shape3))                  stop('Error: Push-shape parameter is not numeric')
  if (!is.numeric(push))                    stop('Error: Push parameter is not numeric')
  if (!is.logical(left))                    stop('Error: Direction input is not a logical value')
  if (!is.logical(include.sd))              stop('Error: include.sd option is not a logical value')
  
  #Check that parameters and options are atomic
  if (length(shape1) != 1)                  stop('Error: Shape1 parameter should be a single number')
  if (length(shape2) != 1)                  stop('Error: Shape2 parameter should be a single number')
  if (length(shape3) != 1)                  stop('Error: Shape3 (push-shape) parameter should be a single number')
  if (length(push)   != 1)                  stop('Error: Push parameter should be a single number')
  if (length(left)   != 1)                  stop('Error: Direction input should be a single logical value')
  if (length(include.sd) != 1)              stop('Error: include.sd option should be a single logical value')
  
  #Check that parameters are in allowable range
  if (min(shape1) <= 0)                     stop('Error: Shape1 parameter should be positive')
  if (min(shape2) <= 0)                     stop('Error: Shape2 parameter should be positive')
  if (min(shape3) < 0)                      stop('Error: Shape3 (push-shape) parameter should be non-negative')
  if (min(push) < 0)                        stop('Error: Push parameter should be between zero and one')
  if (max(push) > 1)                        stop('Error: Push parameter should be between zero and one')
  
  ######################################################################################################
  
  #Compute moments for the left-pushed beta distribution
  if (left) {
    
    #Compute scaling constant (using integration)
    KERNEL <- function(xx) { (xx^(shape1-1))*((1-xx)^(shape2-1))*((1-xx*push)^shape3) }
    INT    <- integrate(f = KERNEL, lower = 0, upper = 1)
    if (INT$message != "OK")   stop('Error: There was an error conducting numerical integration of the density kernel')
    if (INT$value <= 0)        stop('Error: There was an error conducting numerical integration of the density kernel')
    SCALE <- INT$value
    
    #Compute the raw moments
    KERNEL1 <- function(xx) { (xx^(shape1))*((1-xx)^(shape2-1))*((1-xx*push)^shape3) }
    KERNEL2 <- function(xx) { (xx^(shape1+1))*((1-xx)^(shape2-1))*((1-xx*push)^shape3) }
    KERNEL3 <- function(xx) { (xx^(shape1+2))*((1-xx)^(shape2-1))*((1-xx*push)^shape3) }
    KERNEL4 <- function(xx) { (xx^(shape1+3))*((1-xx)^(shape2-1))*((1-xx*push)^shape3) }
    INT1    <- integrate(f = KERNEL1, lower = 0, upper = 1)
    INT2    <- integrate(f = KERNEL2, lower = 0, upper = 1)
    INT3    <- integrate(f = KERNEL3, lower = 0, upper = 1)
    INT4    <- integrate(f = KERNEL4, lower = 0, upper = 1)
    RAW1    <- INT1$value/SCALE
    RAW2    <- INT2$value/SCALE
    RAW3    <- INT3$value/SCALE
    RAW4    <- INT4$value/SCALE
    
    #Compute the central moments
    CENT2 <- RAW2 - RAW1^2
    CENT3 <- RAW3 - 3*RAW1*RAW2 + 2*RAW1^3
    CENT4 <- RAW4 - 4*RAW1*RAW3 + 6*(RAW1^2)*RAW2 - 3*RAW1^4
    
    #Compute the moments
    MEAN <- RAW1
    VAR  <- CENT2
    SKEW <- CENT3/(CENT2^(3/2))
    KURT <- CENT4/(CENT2^2) }
  
  #Compute moments for the left-pushed beta distribution
  if (!left) {
    
    #Compute scaling constant (using integration)
    KERNEL <- function(xx) { (xx^(shape1-1))*((1-xx)^(shape2-1))*((1-push+xx*push)^shape3) }
    INT    <- integrate(f = KERNEL, lower = 0, upper = 1)
    if (INT$message != "OK")   stop('Error: There was an error conducting numerical integration of the density kernel')
    if (INT$value <= 0)        stop('Error: There was an error conducting numerical integration of the density kernel')
    SCALE <- INT$value
    
    #Compute the raw moments
    KERNEL1 <- function(xx) { (xx^(shape1))*((1-xx)^(shape2-1))*((1-push+xx*push)^shape3) }
    KERNEL2 <- function(xx) { (xx^(shape1+1))*((1-xx)^(shape2-1))*((1-push+xx*push)^shape3) }
    KERNEL3 <- function(xx) { (xx^(shape1+2))*((1-xx)^(shape2-1))*((1-push+xx*push)^shape3) }
    KERNEL4 <- function(xx) { (xx^(shape1+3))*((1-xx)^(shape2-1))*((1-push+xx*push)^shape3) }
    INT1    <- integrate(f = KERNEL1, lower = 0, upper = 1)
    INT2    <- integrate(f = KERNEL2, lower = 0, upper = 1)
    INT3    <- integrate(f = KERNEL3, lower = 0, upper = 1)
    INT4    <- integrate(f = KERNEL4, lower = 0, upper = 1)
    RAW1    <- INT1$value/SCALE
    RAW2    <- INT2$value/SCALE
    RAW3    <- INT3$value/SCALE
    RAW4    <- INT4$value/SCALE
    
    #Compute the central moments
    CENT2 <- RAW2 - RAW1^2
    CENT3 <- RAW3 - 3*RAW1*RAW2 + 2*RAW1^3
    CENT4 <- RAW4 - 4*RAW1*RAW3 + 6*(RAW1^2)*RAW2 - 3*RAW1^4
    
    #Compute the moments
    MEAN <- RAW1
    VAR  <- CENT2
    SKEW <- CENT3/(CENT2^(3/2))
    KURT <- CENT4/(CENT2^2) }
  
  #Generate output
  if (include.sd) {
    OUT <- data.frame(mean = 0, var = 0, sd = 0, skew = 0, kurt = 0, excess.kurt = 0) } else {
      OUT <- data.frame(mean = 0, var = 0, skew = 0, kurt = 0, excess.kurt = 0) }
  OUT$mean <- MEAN
  OUT$var  <- VAR
  if (include.sd) { OUT$sd <- sqrt(VAR) }
  OUT$skew <- SKEW
  OUT$kurt <- KURT
  OUT$excess.kurt <- KURT - 3
  
  #Return the output
  OUT }
