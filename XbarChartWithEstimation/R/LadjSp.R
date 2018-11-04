#' LadjSp
#'
#' @examples
#'
#' LadjSp(370,25,5)
#'
#' @import cubature
#'


LadjSp <- function(ARL0nom,m,n) {
  library(cubature)

  if (ARL0nom > 1000 || ARL0nom < 2 || m < 17 || n < 3 || m%%1 != 0 || n%%1 != 0 ) {
    print(paste("Please, revise your entries according to the following conditions:"))
    print(paste("The nominal in-control ARL must be between 2 and 1000"))
    print(paste("The number (m) of Phase I Samples must be equal or larger than 17 and an integer number"))
    print(paste("The size (n) of each Phase I Samples must be equal or larger than 3 and an integer number"))
  }
  else {
    alpha <- 1/ARL0nom
    Lnom <- (-1*qnorm(alpha/2))
    secantc <- function(fun, x0, x1, tol=1e-6, niter=100000){
      for ( i in 1:niter ) {
        funx1 <- fun(x1)
        funx0 <- fun(x0)
        x2 <- ( (x0*funx1) - (x1*funx0) )/( funx1 - funx0 )
        funx2 <- fun(x2)
        if (abs(funx2) < tol) {
          return(x2)
        }
        if (funx2 < 0)
          x0 <- x2
        else
          x1 <- x2
      }
      stop("exceeded allowed number of iteractions")
    }

    ARL0XbarSp <- function (L,m,n) {

      CARL <- function (U) {
        a <- 1/(1 - pnorm(((1/sqrt(m))*qnorm(U[1],0,1))+(L*sqrt(qchisq(U[2],m*(n-1))/(m*(n-1)))),0,1) + pnorm(((1/sqrt(m))*qnorm(U[1],0,1))-(L*sqrt(qchisq(U[2],m*(n-1))/(m*(n-1)))),0,1))
        return(a)
      }
      a <- adaptIntegrate(CARL, lowerLimit = c(0, 0), upperLimit = c(1, 1))$integral
      return (a)
    }

    ARLfunc <- function (alphaf) {
      a <- ARL0XbarSp(alphaf,m,n) - ARL0nom
      return(a)
    }

    cat("This may take several minutes. Please, wait... ")
    L <- secantc(ARLfunc,Lnom*0.9,Lnom*1.1)

    b <- round(ARL0XbarSp(L,m,n),3)

    Lround <- round(L,3)

    print(paste("End of calculations. See results below:"))
    print(paste("L = ", Lround))
    print(paste("When L = ", Lround, ", m = ", m, "and n = ", n, ", the ARL0 = ", b, ", as specified" ))
    print(paste("In summary, this function returned the Limit Factor (L) that generates an Unconditional In-control Average Run Length (ARL0) equal to the specified one for the given number (m) and size (n) of Phase I samples for the Xbar chart with Sp estimator."))
    invisible(L)
  }
}
