#' CDFCARL0Sp
#'
#' @examples
#'
#' CDFCARL0Sp(370,3,25,5)
#'
#' @import cubature
#'


CDFCARL0Sp <- function (t,L,m,n) {
  if (L<=0 ||  L> 5 || m < 1 || n < 2 || m%%1 != 0 || n%%1 != 0 ){
    print(paste("Please, revise your entries according to the following conditions:"))
    print(paste("The Limit Factor should be a posite value equal or smaller than 5"))
    print(paste("The number (m) of Phase I Samples must be equal or larger than 1 and a integer value"))
    print(paste("The size (n) of each Phase I Samples must be equal or larger than 2 and a integer value"))
  }
  else {

    if (t<1){
      print(paste("P( CARL0 <=", t, ")", "= 0 "))
      print(paste("Note that CAR0 is a continuos random variable where CARL0 >=1"))
      print(paste("In Summary, this function returned the Probability of the In-Control Conditional Average Run Length (CARL0) be smaller or equal to", t, "of the specified", L,"-Sigma limits of the Xbar chart for the given number (m) and size (n) of Phase I samples with Sp estimator. This probability is zero, since CAR0 is a continuos random variable where CARL0 >=1"))
      invisible(0)
    }

    else{


      CARL0 <- function (U) {
        a<-pchisq((m*(n-1)*qchisq(1-(1/t), df=1, ncp = (qnorm(U)^2)/m))/(L^2),m*(n-1))
        return(a)
      }

      d <- integrate(CARL0,0,1)$val
      dround <- round(d,4)
      print(paste("P( CARL0 <=", t, ")", "=", dround))
      print(paste("When the Limit Factor (L) = ", L, ", m = ", m, "and n = ", n, ", as specified,", "P(CARL0 <=", t, ")", "=", dround ))
      print(paste("In Summary, this function returned the Probability of the In-Control Conditional Average Run Length (CARL0) be smaller or equal to", t, "of the specified", L,"-Sigma limits of the Xbar chart for the given number (m) and size (n) of Phase I samples with Sp estimator"))
      invisible(d)
    }
  }
}
