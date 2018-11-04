#' ARL0Sp
#'
#' @examples
#'
#' ARL0Sp(3,25,5)
#'
#' @import cubature
#'


ARL0Sp <- function (L,m,n) {
  library(cubature)

  if (L<=0 ||  L> 5 || m < 1 || n < 2 || m%%1 != 0 || n%%1 != 0 ){
    print(paste("Please, revise your entries according to the following conditions:"))
    print(paste("The Limit Factor should be a posite value equal or smaller than 5"))
    print(paste("The number (m) of Phase I Samples must be equal or larger than 1 and a integer value"))
    print(paste("The size (n) of each Phase I Samples must be equal or larger than 2 and a integer value"))
  }
  else {

    CARL <- function (U) {
      a <- 1/(1 - pnorm(((1/sqrt(m))*qnorm(U[1],0,1))+(L*sqrt(qchisq(U[2],m*(n-1))/(m*(n-1)))),0,1) + pnorm(((1/sqrt(m))*qnorm(U[1],0,1))-(L*sqrt(qchisq(U[2],m*(n-1))/(m*(n-1)))),0,1))
      return(a)
    }
    a <- adaptIntegrate(CARL, lowerLimit = c(0, 0), upperLimit = c(1, 1))$integral
    around <- round(a,2)
    print(paste("ARL0 = ", around))
    print(paste("When the Limit Factor (L) = ", L, ", m = ", m, "and n = ", n, ", as specified", ", the ARL0 = ", around ))
    print(paste("In Summary, this function returned the In-Control Unconditional Average Run Length (ARL0) of the specified", L,"-Sigma limits of the Xbar chart for the given number (m) and size (n) of Phase I samples with Sp estimator"))
    invisible(a)
  }
}

