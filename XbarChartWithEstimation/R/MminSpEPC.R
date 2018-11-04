#' MminSpEPC
#'
#' @examples
#'
#' MminSpEPC(370,0.05,0.1,5)
#'
#' @import cubature
#'




MminSpEPC <- function(ARL0nom,e,p,n) {
  library(cubature)

  if ( p > 0.5 || p < 0.01 || e > 0.5 ||e < 0.01 || ARL0nom > 1000 ||ARL0nom < 2 || n < 2  || n%%1 != 0  ) {
    print(paste("Please, revise your entries according to the following conditions:"))
    print(paste("The nominal in-control ARL must be between 2 and 1000"))
    print(paste("The tolerance factor (e) must be between 0.05 and 0.5"))
    print(paste("The probability (p) must be between 0.01 and 0.5"))
    print(paste("The size (n) of each Phase I Samples must be equal or larger than 2 and an integer number"))

  }
  else {
    alpha <- 1/ARL0nom
    L <- (-1*qnorm(alpha/2))


    secant <- function(fun, x0, x1, tol=1e-6, niter=100000){
      for ( i in 1:niter ) {
        funx1 <- fun(x1)
        funx0 <- fun(x0)
        x2 <- ( (x0*funx1) - (x1*funx0) )/( funx1 - funx0 )
        funx2 <- fun(x2)
        if (abs(funx2) < tol) {
          return(x2)
        }
        if (funx2 < 0)
          x1 <- x2
        else
          x0 <- x2
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

    CDFm <- function (m) {
      a <- CDFCARL0((1-e)*ARL0nom,m,n,L) - p
      return(a)
    }

    cat("This may take several minutes. Please, wait... ")
    Tm <- ceiling(secant(CDFm,15,50000))


    r <- round(CDFCARL0((1-e)*ARL0nom,Tm,n,L),2)

    Lr <- round(L,2)

    print(paste("End of calculations. See results below:"))
    print(paste("m = ", Tm))
    print(paste("When L = ", Lr, ", m = ", Tm, "and n = ", n, ", the P(CARL0 <=", (1-e)*ARL0nom, ") =", r,". Note that",(1-e)*ARL0nom,"is e =",e*100,"% smaller than the nominal ARL0, which is",ARL0nom))
    print(paste("In summary, this function returned the minimium number of Phase 1 samples (m) that generates P(CARL0 <=", (1-e)*ARL0nom, ") not larger than", p,"for the",Lr,"- Sigma Xbar chart"," for a given size (n) of the Phase I samples for the Xbar chart with Sp estimator."))
    invisible(Tm)
  }
}
