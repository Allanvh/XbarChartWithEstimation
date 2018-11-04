#' LadjSpEPC
#'
#' @examples
#'
#' LadjSpEPC(370,0.05,0.1,25,5)
#'
#' @import cubature
#'




LadjSpEPC <- function(ARL0nom,e,p,m,n) {
  library(cubature)

  if ( p > 0.5 || p < 0.01 || e > 0.5 ||e < 0 || ARL0nom > 1000 ||ARL0nom < 2 || m < 1 || n < 2 || m%%1 != 0 || n%%1 != 0 ) {
    print(paste("Please, revise your entries according to the following conditions:"))
    print(paste("The nominal in-control ARL must be between 2 and 1000"))
    print(paste("The tolerance factor (e) must be between 0 and 0.5"))
    print(paste("The probability (p) must be between 0.01 and 0.5"))
    print(paste("The number (m) of Phase I Samples must be equal or larger than 1 and an integer number"))
    print(paste("The size (n) of each Phase I Samples must be equal or larger than 2 and an integer number"))

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

    CDFCARL0 <- function (t,m,n,L) {
      CARL0 <- function (U) {
        a<-pchisq((m*(n-1)*qchisq(1-(1/t), df=1, ncp = (qnorm(U)^2)/m))/(L^2),m*(n-1))
        return(a)
      }
      d <- integrate(CARL0,0,1)$val
      return(d)
    }

    CDFaux <- function (s) {
      a <- CDFCARL0((1-e)*ARL0nom,m,n,s)-p
      return (a)
    }
    cat("This may take several minutes. Please, wait... ")
    Lj<-secantc(CDFaux,Lnom*0.9,Lnom*1.1)

    Lround <- round(Lj,3)

    r <- round(CDFCARL0((1-e)*ARL0nom,m,n,Lj),2)

    print(paste("End of calculations. See results below:"))
    print(paste("L = ", Lround))
    print(paste("When L = ", Lround, ", m = ", m, "and n = ", n, ", the P(CARL0 <=", (1-e)*ARL0nom, ") =", r,". Note that",(1-e)*ARL0nom,"is e =",e*100,"% smaller than the nominal ARL0, which is",ARL0nom))
    print(paste("In summary, this function returned the Limit Factor (L) that generates P(CARL0 <=", (1-e)*ARL0nom, ") =", p,"for the given number (m) and size (n) of Phase I samples for the Xbar chart with Sp estimator."))
    invisible(Lj)
  }
}
