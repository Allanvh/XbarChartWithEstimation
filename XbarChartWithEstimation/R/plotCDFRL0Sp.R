#' plotCDFRL0Sp
#'
#' @examples
#'
#' plotCDFRL0Sp(3,25,5)
#'
#' @import cubature
#'



plotCDFRL0Sp <- function (L,m,n) {
  library(cubature)
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

    CDFRL0 <- function (a,m,n,L) {
      t <- floor(a)
      inside <- function (U) {
        a <- 1 - (((pnorm(((1/sqrt(m))*qnorm(U[1],0,1))+(L*sqrt(qchisq(U[2],m*(n-1))/(m*(n-1)))),0,1) - pnorm(((1/sqrt(m))*qnorm(U[1],0,1))-(L*sqrt(qchisq(U[2],m*(n-1))/(m*(n-1)))),0,1))^(t)))
        return(a)
      }
      b <- adaptIntegrate(inside, lowerLimit = c(0, 0), upperLimit = c(1, 1))$integral
      return (b)
    }

    quantileRL0 <-function (p,m,n,L) {
      CDFRL0 <- function (a,m,n,L) {
        t <- a
        inside <- function (U) {
          a <- 1 - (((pnorm(((1/sqrt(m))*qnorm(U[1],0,1))+(L*sqrt(qchisq(U[2],m*(n-1))/(m*(n-1)))),0,1) - pnorm(((1/sqrt(m))*qnorm(U[1],0,1))-(L*sqrt(qchisq(U[2],m*(n-1))/(m*(n-1)))),0,1))^(t)))
          return(a)
        }
        b <- adaptIntegrate(inside, lowerLimit = c(0, 0), upperLimit = c(1, 1))$integral
        return (b)
      }
      CDFm <- function (a) {
        a <- CDFRL0(a,m,n,L) - p
        return(a)
      }
      g<-ceiling(secantc(CDFm,2,1000))
      return(g)
    }
    cat("This may take several minutes. Please, wait... ")
    q <- quantileRL0(0.95,m,n,L)
    q2 <- quantileRL0(0.8,m,n,L)

    par(mar=c(5,5,4,4)+.1)
    CDFRL02 <- Vectorize(CDFRL0)
    curve(CDFRL02(x,m,n,L),0,q ,n=1000,ylim=c(0,1),xlim=c(0,q),xlab="t",ylab="",cex.axis=1.5,type="l",lty=1,yaxs="i",xaxs="i",las=1, lwd=3)
    title(main=paste("P(RL0 <= t)","for", "L=",L, "m=",m, "n=",n ), line=+2.5)

    ARL0nom <- 1/(2*(1-pnorm(L)))
    ARL0nomr <-round(ARL0nom,2)
    CDFARL0nom <- CDFRL0(ARL0nom,m,n,L)
    CDFARL0nomr <- round(CDFARL0nom,2)
    axis(3,ARL0nomr,cex.axis=1,las=1)
    abline(v=ARL0nomr,lty=5.5)
    axis(4,CDFARL0nomr,cex.axis=1,las=1)
    abline(h=CDFARL0nomr,lty=5.5)

    ARL0c <- ARL0(m,n,L)
    CDFmeanrc <- CDFRL0(ARL0c,m,n,L)
    ARL0r <- round(ARL0c,2)
    CDFmeanr <-round(CDFmeanrc,2)
    axis(1,ARL0r,cex.axis=1,las=1, line=1)
    axis(2,CDFmeanr,cex.axis=1,las=1)
    abline(v=ARL0c,lty=5.5,col="blue")
    abline(h=CDFmeanr,lty=5.5,col="blue")

    Median0 <- quantileRL0(0.5,m,n,L)
    CDFmedianc <- CDFRL0(Median0,m,n,L)
    CDFmedianr <-round(CDFmedianc,2)
    axis(1,Median0,cex.axis=1,las=1, line=1)
    axis(2,CDFmedianr,cex.axis=1,las=1)
    abline(v=Median0 ,lty=5.5,col="red")
    abline(h=CDFmedianr,lty=5.5,col="red")
    legend(q2, 0.3, c( paste("MRL0 =", Median0) , paste("ARL0 =", ARL0r), paste("Nominal ARL0 =", ARL0nomr)), cex=1, lty=c(5.5,5.5,5.5),lwd=c(1,1,1),col=c("red","blue","black"))
    print(paste("End of calculations."))
    print(paste("This function returned the plot of the CDF of the In-Control Unconditional Run Length (RL0)."))
  }
}

