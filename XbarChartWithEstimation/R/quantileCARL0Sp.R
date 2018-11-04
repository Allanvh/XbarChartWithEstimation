#' quantileCARL0Sp
#'
#' @examples
#'
#' quantileCARL0Sp(0.5,3,25,5)
#'
#' @import cubature
#'




quantileCARL0Sp <-function (p,L,m,n) {

  library(cubature)

  if (p > 0.95 || p < 0.05 ||L<2 || L> 4  || m < 20 || n < 3 || m%%1 != 0 || n%%1 != 0 ){
    print(paste("Please, revise your entries according to the following conditions:"))
    print(paste("The specified quantile must be between 0.05 and 0.95"))
    print(paste("The Limit Factor should be equal or larger than 2 and equal or samller than 4"))
    print(paste("The number (m) of Phase I Samples must be equal or larger than 20 and a integer value"))
    print(paste("The size (n) of each Phase I Samples must be equal or larger than 3 and a integer value"))
  }
  else {

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



    CDFm <- function (a) {
      a <- CDFCARL0(a,m,n,L) - p
      return(a)
    }

    cat("This may take several minutes. Please, wait... ")
    g<-secantc(CDFm,10,1000)

    ground <- round(g,2)

    print(paste("End of calculations. See results below:"))
    print(paste("P( CARL0 <=",ground,") =", p))
    print(paste("The",p,"-quantile of CARL0 is ", ground))
    print(paste("When the Limit Factor (L) = ", L, ", m = ", m, "and n = ", n, ", as specified, the", p , "-quantile of the CARL0 is ", ground ))
    print(paste("In Summary, this function returned the",p,"-quantile of the In-Control Conditional Average Run Length (RL0) of the specified", L,"-Sigma limits of the Xbar chart for the given number (m) and size (n) of Phase I samples with Sp estimator"))
    invisible(g)
  }
}

