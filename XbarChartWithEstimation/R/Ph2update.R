check.packages <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
        install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
}


###########################################################################################

packages<-c("devtools")
install_github('bolus123/PH1AND2XBAR')

###########################################################################################

#require(mvtnorm)

#source(file = 'https://raw.githubusercontent.com/bolus123/PH1AND2XBAR/master/R/main.R')

###########################################################################################

PH2.L <- PH2.get.cc
