#' CDFRL0Sp
#'
#' @examples
#'
#' CDFRL0Sp(370,3,25,5)
#'
#' @import cubature
#'


CDFRL0Sp <- function (a,L,m,n) {
  library(cubature)
  if (L<=0 || L> 5 ||  m < 1 || n < 2 || m%%1 != 0 || n%%1 != 0 ){
    print(paste("Please, revise your entries according to the following conditions:"))
    print(paste("The Limit Factor should be a posite value equal or smaller than 5"))
    print(paste("The number (m) of Phase I Samples must be equal or larger than 1 and a integer value"))
    print(paste("The size (n) of each Phase I Samples must be equal or larger than 2 and a integer value"))
  }

  else {

    if (a<1){
      print(paste("P( RL0 <=", a, ")", "= 0 "))
      print(paste("Note that R0 is a discrete random variable where RL0 >=1"))
      print(paste("In Summary, this function returned the Probability of the In-Control Unconditional Run Length (RL0) be smaller or equal to", a, "of the specified", L,"-Sigma limits of the Xbar chart for the given number (m) and size (n) of Phase I samples with Sp estimator. This probability is zero, since R0 is a discrete random variable where RL0 >=1"))
      invisible(0)
    }

    else {
      t <- floor(a)
      inside <- function (U) {
        a <- 1 - (((pnorm(((1/sqrt(m))*qnorm(U[1],0,1))+(L*sqrt(qchisq(U[2],m*(n-1))/(m*(n-1)))),0,1) - pnorm(((1/sqrt(m))*qnorm(U[1],0,1))-(L*sqrt(qchisq(U[2],m*(n-1))/(m*(n-1)))),0,1))^(t)))
        return(a)
      }
      b <- adaptIntegrate(inside, lowerLimit = c(0, 0), upperLimit = c(1, 1))$integral
      bround <- round(b,4)
      print(paste("P( RL0 <=", a, ")", "=", bround))
      print(paste("When the Limit Factor (L) = ", L, ", m = ", m, "and n = ", n, ", as specified,", "P(RL0 <=", a, ")", "=", bround ))
      print(paste("In Summary, this function returned the Probability of the In-Control Unconditional Run Length (RL0) be smaller or equal to", a, "of the specified", L,"-Sigma limits of the Xbar chart for the given number (m) and size (n) of Phase I samples with Sp estimator"))
      invisible(b)
    }
  }
}

