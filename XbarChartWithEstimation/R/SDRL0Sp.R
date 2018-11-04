#' SDRL0Sp
#'
#' @examples
#'
#' SDRL0Sp(3,25,5)
#'
#' @import cubature
#'


SDRL0Sp <- function (L,m,n) {
  library(cubature)

  if (L<=0 || L> 5 || m < 25 || n < 5 || m%%1 != 0 || n%%1 != 0 ){
    print(paste("Please, revise your entries according to the following conditions:"))
    print(paste("The Limit Factor should be a posite value equal or smaller than 5"))
    print(paste("The number (m) of Phase I Samples must be equal or larger than 25 and a integer value"))
    print(paste("The size (n) of each Phase I Samples must be equal or larger than 5 and a integer value"))
  }
  else {

    ARL0 <- function (m,n,L) {
      CARL <- function (U) {
        a <- 1/(1 - pnorm(((1/sqrt(m))*qnorm(U[1],0,1))+(L*sqrt(qchisq(U[2],m*(n-1))/(m*(n-1)))),0,1) + pnorm(((1/sqrt(m))*qnorm(U[1],0,1))-(L*sqrt(qchisq(U[2],m*(n-1))/(m*(n-1)))),0,1))
        return(a)
      }
      a <- adaptIntegrate(CARL, lowerLimit = c(0, 0), upperLimit = c(1, 1))$integral
      return (a)
    }

    ARL02 <- function (m,n,L) {

      CARL <- function (U) {
        a <- (1/(1 - pnorm(((1/sqrt(m))*qnorm(U[1],0,1))+(L*sqrt(qchisq(U[2],m*(n-1))/(m*(n-1)))),0,1) + pnorm(((1/sqrt(m))*qnorm(U[1],0,1))-(L*sqrt(qchisq(U[2],m*(n-1))/(m*(n-1)))),0,1)))^2
        return(a)
      }
      a <- adaptIntegrate(CARL, lowerLimit = c(0, 0), upperLimit = c(1, 1))$integral
      return (a)

    }

    VCARL0 <- function (m,n,L) {
      a <- ARL02(m,n,L) - (ARL0(m,n,L))^2
      return (a)
    }

    SDCARL <- function (m,n,L) {
      a <- sqrt( ARL02(m,n,L) - (ARL0(m,n,L))^2)
      return (a)
    }

    AVCRL0 <- function (m,n,L) {
      CFAR <- function (U) {
        a <- (1 - pnorm(((1/sqrt(m))*qnorm(U[1],0,1))+(L*sqrt(qchisq(U[2],m*(n-1))/(m*(n-1)))),0,1) + pnorm(((1/sqrt(m))*qnorm(U[1],0,1))-(L*sqrt(qchisq(U[2],m*(n-1))/(m*(n-1)))),0,1))
        b <- (1-a)/(a^2)
        return(b)
      }
      a <- adaptIntegrate(CFAR, lowerLimit = c(0, 0), upperLimit = c(1, 1))$integral
      return (a)
    }

    VRL0 <- function (m,n,L) {
      a <- AVCRL0(m,n,L) + VCARL0(m,n,L)
      return (a)
    }

    a <- sqrt(VRL0(m,n,L))
    around <- round(a,2)
    print(paste("SDRL0 = ", around))
    print(paste("When the Limit Factor (L) = ", L, ", m = ", m, "and n = ", n, ", as specified", ", the SDRL0 = ", around ))
    print(paste("In Summary, this function returned the standar deviation of the In-Control Unconditional Run Length (RL0) of the specified", L,"-Sigma limits of the Xbar chart for the given number (m) and size (n) of Phase I samples with Sp estimator"))
    invisible(a)
  }
}
