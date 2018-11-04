#' plotPDFCARL0
#'
#' @examples
#'
#' plotPDFCARL0(3,25,5)
#'
#' @import cubature
#' @import numDeriv
#'


plotPDFCARL0 <- function (L,m,n) {
  library(cubature)
  library(numDeriv)

  if (L<2 || L> 4 || m < 5 || n < 3 || m%%1 != 0 || n%%1 != 0 ){
    print(paste("Please, revise your entries according to the following conditions:"))
    print(paste("The Limit Factor should be equal or larger than 2 and equal or samller than 4"))
    print(paste("The number (m) of Phase I Samples must be equal or larger than 5 and a integer value"))
    print(paste("The size (n) of each Phase I Samples must be equal or larger than 3 and a integer value"))
  }
  else{

    dev.new()

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

    ARL0 <- function (m,n,L) {
      CARL <- function (U) {
        a <- 1/(1 - pnorm(((1/sqrt(m))*qnorm(U[1],0,1))+(L*sqrt(qchisq(U[2],m*(n-1))/(m*(n-1)))),0,1) + pnorm(((1/sqrt(m))*qnorm(U[1],0,1))-(L*sqrt(qchisq(U[2],m*(n-1))/(m*(n-1)))),0,1))
        return(a)
      }
      a <- adaptIntegrate(CARL, lowerLimit = c(0, 0), upperLimit = c(1, 1))$integral
      return (a)
    }


    CDFCARL0 <- function (t,m,n,L) {
      CARL0 <- function (U) {
        a<-pchisq((m*(n-1)*qchisq(1-(1/t), df=1, ncp = (qnorm(U)^2)/m))/(L^2),m*(n-1))
        return(a)
      }
      d <- integrate(CARL0,0,1)$val
      return(d)
    }

    quantileCARL0 <-function (p,m,n,L) {
      CDFm <- function (a) {
        a <- CDFCARL0(a,m,n,L) - p
        return(a)
      }
      g<-secantc(CDFm,10,1000)
      return(g)
    }

    CDF <- function (h) {
      g <- CDFCARL0(h,m,n,L)
      return(g)
    }
    PDF <- function (x) {
      f <- grad(CDF, x)
      return(f)
    }

    cat("This may take several minutes. Please, wait... ")

    q <- quantileCARL0(0.95,m,n,L)
    q2 <- quantileCARL0(0.8,m,n,L)
    ARL0nom <- 1/(2*(1-pnorm(L)))
    ARL0nomr <-round(ARL0nom,2)

    par(mar=c(5,5,4,4)+.1)

    PDF2 <- Vectorize(PDF)
    curve(PDF2,1.1,q,xlim=c(0,q),xlab="t",ylab="",n=100,cex.axis=1.5,type="l",lty=1,lwd=3,las=1)
    title(main=paste("pdf of the IC CARL","for", "L=",L, "m=",m, "n=",n  ), line=+2.5)
    axis(3,ARL0nomr,cex.axis=1,las=1)
    abline(v=ARL0nomr,lty=5.5)

    ARL0r <- round(ARL0(m,n,L),2)
    axis(1,ARL0r,cex.axis=1,las=1,line=1)
    abline(v=ARL0(m,n,L),lty=5.5,col="blue")

    Median0 <- round(quantileCARL0(0.5,m,n,L),2)
    axis(1,Median0,cex.axis=1,las=1,line=1)
    abline(v=Median0 ,lty=5.5,col="red")

    print(paste("End of calculations."))
    print(paste("This function returned the plot of the PDF of the In-Control Conditional Average Run Length (CARL0)."))
  }
}
