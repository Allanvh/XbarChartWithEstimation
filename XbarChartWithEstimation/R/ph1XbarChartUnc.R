#' ph1XbarChartUnc
#'
#' @examples
#'
#' ph1XbarChartUnc()
#'
#' @import mvtnorm
#'




ph1XbarChartUnc <- function() {

  
  cat("Please save your Excel table with your Phase I data as .txt with only the main data. I.e., with no header, columm or row names. Each of the Phase I samples must be in a diffrent row. The number (m) of Phase I samples shoud be larger than 16 and the size (n) of each sample shoud be larger than 2. Finally, your data should be i.i.d. normally distributed. If your data follows all these conditions, plase write 'yes' and select your .txt file:\n\n")
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
    print(paste("The grand mean estimate (Xbarbar) = ", xbarbar))
    print(paste("The pooled standard deviation estimate (Sp) = ", Sp))
    
    cat("\n")
    cat("Please, choose a target (less than 0.5) for the Joint False Alarm Probability (FAP). A commom values is 0.05. \n\n")
    ARL0nom2 <- readline()
    
    ARL0nom <- as.numeric(ARL0nom2)
    
    if  (ARL0nom > 0.5 ) {
      
      stop("Joint False Alarm Probability (FAP) must be smaller or equal to 0.5")
      
    }
    
    else {
      
      cat("\n")
      cat("The Phase I control limits is being calculated. Please, wait...")
      cat("\n")
      cat("\n")
      
      ph1LadjSp <- function(a,b,c) {
        library(mvtnorm)
        ph1LadjSpt <- function(
          k
          ,nu
          ,FAP = 0.1
          ,off.diag = -1/(k - 1)
          #,alternative = '2-sided'
          ,c4.option = FALSE
          ,maxiter = 10000
          ,method = 'direct'
          ,indirect.interval = c(1, 7)
          ,indirect.subdivisions = 100L
          ,indirect.tol = .Machine$double.eps^0.25
          
          
        ){
          
          ####################################################################################################################################################
          # Phase I Xbar Chart Based on Champ and Jones(2004); see also Yao et al. (2017), and Yao and Chakraborti (2018)
          ####################################################################################################################################################
          
          #use too many functions from package mvtnorm
          require(mvtnorm)
          #require(adehabitatLT)
          #dchi <- adehabitatLT::dchi
          
          #focus on these 3 functions from package mvtnorm
          
          
          #pmvt <- mvtnorm::pmvt
          #qmvt <- mvtnorm::qmvt
          #pmvnorm <- mvtnorm::pmvnorm
          
          
          ####################################################################################################################################################
          
          #source('https://raw.githubusercontent.com/bolus123/Statistical-Process-Control/master/MKLswitch.R')
          
          ####################################################################################################################################################
          #parts of getting the charting constant l
          ####################################################################################################################################################
          
          c4.f <- function(nu) sqrt(2 / nu) * 1 / beta(nu / 2, 1 / 2) * sqrt(pi)             #c4.function
          
          corr.par.f <- function(nu, c4.option) {
            if (c4.option == TRUE) {
              corr.par <- c4.f(nu)
            } else {
              corr.par <- 1
            }
          }
          
          PH1.corr.f <- function(k, off.diag = - 1 / (k - 1)){
            #correlation matrix
            crr <- diag(k)
            crr[which(crr == 0)] <- off.diag
            
            crr
            
          }
          
          #first.der.c4.f <- function(nu){
          #
          #    beta.part <- (beta(nu / 2, 1 / 2)) ^ (-1)
          #
          #    digamma1 <- digamma(nu / 2)
          #    digamma2 <- digamma(nu / 2 + 1 / 2)
          #
          #    - sqrt(2 * pi) / 2 * beta.part * ((digamma1 - digamma2) / sqrt(nu) + nu ^ (- 3 /2))
          #
          #}
          #
          #first.der.cons.f <- function(nu, tau){
          #
          #    1 / 2 * (1 - 1 / (nu - tau)) ^ (- 1 / 2) * (nu - tau) ^ - 2
          #
          #}
          #
          #cons.f <- function(nu, tau){
          #
          #    sqrt((nu - tau - 1) / (nu - tau))
          #
          #}
          
          
          ####################################################################################################################################################
          #get l using the multivariate t cdf
          ####################################################################################################################################################
          
          PH1.get.cc.mvt <- function(
            k
            ,nu
            ,FAP = 0.1
            #,Phase1 = TRUE
            ,off.diag = -1/(k - 1)
            #,alternative = '2-sided'
            ,c4.option = TRUE
            ,maxiter = 10000
            
          ){
            
            alternative = '2-sided'                                                   #turn off the alternative
            #The purpose of this function is
            #MCMC <- FALSE                                                             #to obtain l and k based on
            #multivariate t.
            #MCMC part is not available now.
            
            #if (off.diag == NULL) off.diag <- ifelse(Phase1 == TRUE, - 1 /(m - 1), 1 / (m + 1))
            
            corr.P <- PH1.corr.f(k = k, off.diag = off.diag)                              #get correlation matrix with equal correlations
            
            pu <- 1 - FAP
            
            corr.par <- corr.par.f(nu, c4.option)
            
            #if (MCMC == TRUE) {
            #MVN.Q.Gibbs.Sampling(
            #    pu,
            #    MCMC.maxsim,
            #    corr.P,
            #    tails = alternative,
            #    burn = MCMC.burn,
            #    search.interval = MCMC.search.interval
            #)
            
            
            
            #} else {
            
            L <- ifelse(
              alternative == '2-sided',
              qmvt(pu, df = nu, sigma = corr.P, tail = 'both.tails', maxiter = maxiter)$quantile,
              qmvt(pu, df = nu, sigma = corr.P, maxiter = maxiter)$quantile
            )
            #get L by multivariate T
            
            
            
            #}
            
            
            
            c.i <- L * corr.par * sqrt((k - 1) / k)             #get K
            
            list(l = L, c.i = c.i)
            
          }
          
          ####################################################################################################################################################
          #get L by multivariate Normal
          ####################################################################################################################################################
          
          #Old code is using chi distribution from package adehabitatLT
          #
          #PH1.joint.pdf.mvn.chisq <- function(Y, K, m, nu, sigma, alternative = '2-sided') {
          #
          #    s <- length(Y)
          #
          #    L <- K / sqrt((m - 1) / m * nu) * Y / c4.f(nu)
          #
          #    dpp <- lapply(
          #            1:s,
          #            function(i){
          #
          #                LL <- rep(L[i], m)
          #
          #                ifelse(
          #                    alternative == '2-sided',
          #                    pmvnorm(lower = -LL, upper = LL, sigma = sigma),
          #                    pmvnorm(lower = rep(-Inf, m), upper = LL, sigma = sigma)
          #                )
          #
          #            }
          #
          #    )
          #
          #    dpp <- unlist(dpp)
          #
          #    dpp * dchi(Y, df = nu)
          #
          #
          #}
          
          PH1.joint.pdf.mvn.chisq <- function(
            Y
            ,c.i
            ,k
            ,nu
            ,sigma
            #,alternative = '2-sided'
            ,c4.option = TRUE)
          {
            
            alternative = '2-sided'                                                   #turn off the alternative
            
            s <- length(Y)
            
            corr.par <- corr.par.f(nu, c4.option)
            
            L <- c.i / sqrt((k - 1) / k * nu) * sqrt(Y) / corr.par
            
            dpp <- lapply(
              1:s,
              function(i){
                
                LL <- rep(L[i], k)
                
                ifelse(
                  alternative == '2-sided',
                  pmvnorm(lower = -LL, upper = LL, sigma = sigma),
                  pmvnorm(lower = rep(-Inf, k), upper = LL, sigma = sigma)
                )
                
              }
              
            )
            
            dpp <- unlist(dpp)
            
            dpp * dchisq(Y, df = nu)
            
            
          }
          
          PH1.root.mvn.F <- function(
            c.i
            , k
            , nu
            , sigma
            , pu
            #, alternative = '2-sided'
            , c4.option = TRUE
            , subdivisions = 2000
            , rel.tol = 1e-2)
          {
            alternative = '2-sided'                                                   #turn off the alternative
            #The purpose of this function is
            #s <- length(Y)                                          #to obtain appropriate L
            #                                                        #by multivariate normal
            #L <- K / sqrt((m - 1) / m * nu) * Y * c4.f(nu)          #
            #
            #pp <- lapply(
            #        1:s,
            #        function(i){
            #
            #            LL <- rep(L[i], m)
            #            ifelse(
            #                alternative == '2-sided',
            #                pmvnorm(lower = -LL, upper = LL, sigma = sigma),
            #                pmvnorm(lower = rep(-Inf, m), upper = LL, sigma = sigma)
            #            )
            #
            #        }
            #
            #)
            #
            #pp <- mean(unlist(pp))
            
            
            pp <- integrate(
              PH1.joint.pdf.mvn.chisq,
              lower = 0,
              upper = Inf,
              c.i = c.i,
              k = k,
              nu = nu,
              sigma = sigma,
              #alternative = alternative,
              c4.option = c4.option,
              subdivisions = subdivisions,
              rel.tol = rel.tol
            )$value
            
            pu - pp
            
            
          }
          
          
          PH1.get.cc.mvn <- function(
            k
            ,nu
            ,FAP = 0.1
            #,Phase1 = TRUE
            ,off.diag = -1/(k - 1)
            #,alternative = '2-sided'
            ,c4.option = TRUE
            ,interval = c(1, 7)
            ,maxiter = 10000
            ,subdivisions = 2000
            ,tol = 1e-2
            
          ){
            alternative = '2-sided'                                                   #turn off the alternative
            #The purpose of this function is
            #to obtain L and K based on
            #MCMC <- FALSE                                           #multivariate normal.
            #MCMC part is not available now.
            #if (is.null(off.diag)) off.diag <- ifelse(Phase1 == TRUE, - 1 /(m - 1), 1 / (m + 1))
            
            corr.par <- corr.par.f(nu, c4.option)
            
            corr.P <- PH1.corr.f(k = k, off.diag = off.diag)
            
            pu <- 1 - FAP
            
            #Y <- sqrt(rchisq(maxsim, df = nu))
            
            #if (MCMC == TRUE) {
            
            #MVN.Q.Gibbs.Sampling(
            #    pu,
            #    MCMC.maxsim,
            #    corr.P,
            #    tails = alternative,
            #    burn = MCMC.burn,
            #    search.interval = MCMC.search.interval
            #)
            
            
            
            #} else {
            
            
            c.i <- uniroot(
              PH1.root.mvn.F,
              interval = interval,
              k = k,
              nu = nu,
              #Y = Y,
              sigma = corr.P,
              pu = pu,
              #alternative = alternative,
              c4.option = c4.option,
              subdivisions = subdivisions,
              tol = tol,
              rel.tol = tol,
              maxiter = maxiter
            )$root
            
            #}
            
            L <- c.i / corr.par * sqrt(k / (k - 1))
            
            list(l = L, c.i = c.i)
            
            
          }
          
          #get.cc.mvn(20, 80)
          
          ####################################################################################################################################################
          
          
          
          alternative = '2-sided'                                                   #turn off the alternative
          #The purpose of this function is to obtain L and K
          #by multivariate T or multivariate normal.
          #Multivariate normal is so time-consuming
          #that I do not recommend.
          
          #Phase1 <- TRUE                                  #need to delete if need to use Phase 2
          
          #method <- 'direct'								#need to delete if need to use indirect
          
          #if (is.null(off.diag)) off.diag <- ifelse(Phase1 == TRUE, - 1 /(m - 1), 1 / (m + 1))
          
          is.int <- ifelse(nu == round(nu), 1, 0)
          
          #if (is.int == 0) cat('Nu is not an integrer. Please check.', '\n') #need to delete if need to use indirect
          
          #if (method == 'direct' & is.int == 1) {                       #using multivariate T to obtain L and K
          
          if (is.int == 1 && method == 'direct') {
            #need to delete if need to use indirect
            PH1.get.cc.mvt(
              k = k
              ,nu = nu
              ,FAP = FAP
              #,Phase1 = Phase1
              ,off.diag = off.diag
              #,alternative = alternative
              ,c4.option = c4.option
              ,maxiter = maxiter
            )
            
          } else if (is.int == 1 && method == 'indirect'){
            
            cat('Nu is an integer. Using the indirect method may slow the computation process down.', '\n')
            
            PH1.get.cc.mvn(
              k = k
              ,nu = nu
              ,FAP = FAP
              #,Phase1 = Phase1
              ,off.diag = off.diag
              #,alternative = alternative
              ,c4.option = c4.option
              ,interval = indirect.interval
              #,maxsim = indirect.maxsim
              ,subdivisions = indirect.subdivisions
              ,maxiter = maxiter
              ,tol = indirect.tol
            )
            
          }
          
          else if (is.int == 0 && method == 'direct'){   #using multivariate normal to obtain L and K
            
            
            stop('Nu is not an integer. Please use the indirect method instead.')
            
            
          } else if (is.int == 0 && method == 'indirect') {
            
            PH1.get.cc.mvn(
              k = k
              ,nu = nu
              ,FAP = FAP
              #,Phase1 = Phase1
              ,off.diag = off.diag
              #,alternative = alternative
              ,c4.option = c4.option
              ,interval = indirect.interval
              #,maxsim = indirect.maxsim
              ,subdivisions = indirect.subdivisions
              ,maxiter = maxiter
              ,tol = indirect.tol
            )
            
          } else {
            
            stop('Unknown method. Please select the direct method or the indirect method.')
            
          }
          
          
        }
        
        d <- ph1LadjSpt(a,a*(b-1),c)$c.i
        return (d)
        
      }
      
      Ladj<- ph1LadjSp(m,n,ARL0nom)
      UCL <- xbarbar + Ladj*(Sp/sqrt(n))
      LCL <- xbarbar - Ladj*(Sp/sqrt(n))
      
      UCL2 <- round(UCL,2)
      LCL2 <- round(LCL,2)
      CL2 <- round(xbarbar,2) 
      
      print(paste("The Phase I Upper Control Limit (UCL) = ", UCL2))
      print(paste("The Central Line (CL) = ", CL2))
      print(paste("The Phase I Lower Control Limit (LCL) = ", LCL2))
      
  
      
      dev.new()
      par(mar=c(5,5,4,4)+.1)
      
      
      t <- seq(from = 1, to = m, by = 1)
    
      plot(t,xbar,xlim=c(1,m),ylim=c(pmin(min(xbar),LCL2),pmax(max(xbar),UCL2)),xlab="Xbar",ylab="",cex.axis=1.5,lty=2,lwd=5,yaxt='n',xaxt='n',pch=16)
      axis(1,at=t,labels=t,cex.axis=1.5,las=1)
      
      
      
      abline(h=UCL2,lty=1,,lwd=3)
      axis(2,UCL2,cex.axis=1,las=1)
      abline(h=LCL2,lty=1,,lwd=3)
      axis(2,LCL2,cex.axis=1,las=1)
      abline(h=CL2,lty=1,,lwd=3)
      axis(2,CL2,cex.axis=1,las=1)
      ARL0nom2 <- CL2 <- round(ARL0nom,2) 
      
      if ((min(xbar)<LCL2)||(max(xbar)>UCL2)) {
        
        title(main=paste("Phase I Xbar Control Chart","for", "FAP = ",ARL0nom2, "m = ",m, "n = ",n , "\nWarning: The data may came from an out-of-control process!"))
        
      }
       else {
         title(main=paste("Phase I Xbar Control Chart","for", "FAP = ",ARL0nom2, "m = ",m, "n = ",n))
         
         
       }
     
    }
  }
  
  else {
    stop("Please, do a Phase I analysis. Your data must satisfies all the assumptions!")
  }
  
}


