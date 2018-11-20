XbarChartUnc <- function() {
  library(cubature)
  
  cat("Please save your Excel table as a .txt with only the main data. I.e., with no header, columm and row names. Each of the Phase I samples must be in a diffrent row. The number (m) of Phase I samples shoud be larger than 16 and the size (n) of each sample shoud be larger than 2. Finally, your data should came fom an in-control process and be i.i.d. normally distributed. If your data follows all of these assumpitions, plase write 'yes' and select your .txt file:\n\n")
  y <- readline()
  
  if (y == "yes" | y == "'yes'" | y == "'y'" | y == "y") {
    
    dev.new()
    data <- read.table(file.choose())
    X <- data.matrix(data)
    
    m <- nrow(X)
    n <- ncol(X)
    
    if (m < 17 || n < 3 ) {
      stop("Impossible to continue because the number (m) of Phase I samples shoud be larger than 16 and the size (n) of each sample shoud be larger than 2. Please get more Phase I data and came back.")
    }
    
    xbar <- rep(NA, m)
    for (i in 1:m){     
      xbar[i] = mean(X[i,])
    }
    xbarbar <- mean(xbar)
    
    S2 <- rep(NA, m)
    for (i in 1:m){     
      S2[i] = var(X[i,])
    }
    S2p <- mean(S2)
    Sp <- sqrt(S2p)
    
    print(paste("Number of Phase I samples (m) = ", m))
    print(paste("Size of each Phase I samples (n) = ", n))
    print(paste("The grand mean estimate (Xbarbar) = ", round(xbarbar,2)))
    print(paste("The pooled standard deviation estimate (Sp) = ", round(Sp,2)))
    
    cat("\n")
    cat("Please, choose a target (between 2 and 1000) for the In-control Unconditional Average Run Lenght (ARL0). The most commom values is 370. \n\n")
    ARL0nom2 <- readline()
    
    ARL0nom <- as.numeric(ARL0nom2)
    
    if  (ARL0nom > 1000 || ARL0nom < 2 ) {
      
      stop("The In-control Unconditional Average Run Lenght (ARL0) must be between 2 and 1000.")
      
    }
    
    else {
      
      cat("\n")
      cat("The unconditional adjusted control limits is being calculated. Please, wait...")
      cat("\n")
      cat("\n")
      
      LadjSp <- function(ARL0nom,m,n) {
        
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
        
        L <- secantc(ARLfunc,Lnom*0.9,Lnom*1.1)
        
        b <- round(ARL0XbarSp(L,m,n),3)
        
        Lround <- round(L,3)
        
        return(L)
      }
      
      Ladj<- LadjSp(ARL0nom,m,n)
      UCL <- xbarbar + Ladj*(Sp/sqrt(n))
      LCL <- xbarbar - Ladj*(Sp/sqrt(n))
      
      UCL2 <- round(UCL,2)
      LCL2 <- round(LCL,2)
      CL2 <- round(xbarbar,2) 
      
      print(paste("The Unconditional Upper Control Limit (UCL) = ", UCL2))
      print(paste("The Central Line (CL) = ", CL2))
      print(paste("The Unconditional Lower Control Limit (LCL) = ", LCL2))
      
      
      
      dev.new()
      par(mar=c(5,5,4,4)+.1)
      plot.new(); plot.window(xlim=c(1,10),ylim=c(pmin(min(xbar),LCL2),pmax(max(xbar),UCL2)) )
      abline(h=UCL2,lty=1,lwd=3)
      axis(2,UCL2,cex.axis=1,las=1)
      abline(h=LCL2,lty=1,lwd=3)
      axis(2,LCL2,cex.axis=1,las=1)
      abline(h=CL2,lty=1,lwd=3)
      axis(2,CL2,cex.axis=1,las=1)
      ARL0nom2 <- CL2 <- round(ARL0nom,2) 
      title(main=paste("Xbar Control Chart","for", "ARL0 = ",ARL0nom2, "m = ",m, "n = ",n ))
      
      cat("\n")
      cat("\n")
      cat("Do you have a Phase II Sample? If so, write Yes\n\n")
      
      y <- readline()
      
      if (y == "yes" | y == "'yes'" | y == "'y'" | y == "y") {
        
        x2 <- rep(NA, n)
        
        for(j in 1:n) {
          
          cat("\n")
          cat("Enter your Phase II data number",j,":\n\n")
          gh <- readline()
          x2[j] <- as.numeric(gh)
        
        }
        
        h <- mean(x2)
        
        points(1, h)
      
      }
    
      else{
        
        cat("\n")
        cat("End of the Program.")
        
      }
      
      }
  }
  
  else {
    stop("Please, do a Phase I analysis. Your data must satisfies all the assumptions!")
  }
  
}

