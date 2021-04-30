####### This R code will calculate N_DL and N_MN for the worked example from the multinomial sample size paper


#############################################################################################################
##### Calculate N_DL which assumes an R2_NAGELKERKE of 0.15 for each distinct logistic regression model #####
#############################################################################################################

### First define the number of events in each category
EV1 <- 800
EV2 <- 55
EV3 <- 169
EV4 <- 42

### Next write a function to calculate Lnull for a logistic regression model with two outcome category prevalences
### Also output the number of events involved in calculating the lnull, as this is also required to calculating max(R2_CS)
Lnull.calc.DL <- function(E1, E2){output.lnull <- E1*log(E1/(E1 + E2)) + (E2)*log(E2/(E1+E2))
                               output.n <- E1 + E2
                               output <- c(output.lnull, output.n)
                               names(output) <- c("lnull","n")
                               return(output)
}

### Calculate lnull for each distinct logistic regression model
Lnull2 <- Lnull.calc.DL(EV1,EV2)
Lnull3 <- Lnull.calc.DL(EV1,EV3)
Lnull4 <- Lnull.calc.DL(EV1,EV4)

### Calculate max(R2_CS)
maxR2.2 <- 1 - exp(2*Lnull2["lnull"]/Lnull2["n"])
maxR2.3 <- 1 - exp(2*Lnull3["lnull"]/Lnull3["n"])
maxR2.4 <- 1 - exp(2*Lnull4["lnull"]/Lnull4["n"])

maxR2.2
maxR2.3
maxR2.4

### Calculate R2_CS which assumes an R2_NAGELKKERKE of 0.15
R2_CS_NAGEL.2 <- 0.15*maxR2.2
R2_CS_NAGEL.3 <- 0.15*maxR2.3
R2_CS_NAGEL.4 <- 0.15*maxR2.4

### Calculate N_DL_NAGELKERKE based on this assumed R2_CS value, to target a threshold of 0.9, by plugging into formula of Riley et al.
## Let p be number of predictors considered in initial model, supposing we will consider the same set before variable selection
p <- 16

## Let S be the level of shrinkage we are targeting
S <- 0.9

### Calculate N_DL_NAGELKERKE temporary
N_DL_NAGEL.2.temp <- p/((S - 1)*log(1 - R2_CS_NAGEL.2/S))
N_DL_NAGEL.3.temp <- p/((S - 1)*log(1 - R2_CS_NAGEL.3/S))
N_DL_NAGEL.4.temp <- p/((S - 1)*log(1 - R2_CS_NAGEL.4/S))

N_DL_NAGEL.2.temp
N_DL_NAGEL.3.temp
N_DL_NAGEL.4.temp

### Note that this is the number of individuals required for each distinct logistic regression model, but only a subset of the individuals
### recruited will be used in these models, so need to divide by the proportion of individuals that are used in each of the models
n.total <- EV1 + EV2 + EV3 + EV4
Eprop.total.2 <- (EV1 + EV2)/n.total
Eprop.total.3 <- (EV1 + EV3)/n.total
Eprop.total.4 <- (EV1 + EV4)/n.total

### Calculate N_DL_NAGELKERKE final
N_DL_NAGEL.2 <- N_DL_NAGEL.2.temp/Eprop.total.2
N_DL_NAGEL.3 <- N_DL_NAGEL.3.temp/Eprop.total.3
N_DL_NAGEL.4 <- N_DL_NAGEL.4.temp/Eprop.total.4

N_DL_NAGEL.2.temp
N_DL_NAGEL.3.temp
N_DL_NAGEL.4.temp

N_DL_NAGEL.2
N_DL_NAGEL.3
N_DL_NAGEL.4


###################################################################################################
##### Calculate N_MN which assumes an R2_NAGELKERKE of 0.15 for the overall multinomial model #####
###################################################################################################

### First write a function to calculate Lnull for a multinomial model
### Also output the number of events involved in calculating the lnull, as this is also required to calculating max(R2_CS)
Lnull.calc.MN <- function(E1, E2, E3, E4){n <- E1 + E2 + E3 + E4
                                  output.lnull <- E1*log(E1/n) + (E2)*log(E2/n) + (E3)*log(E3/n) + (E4)*log(E4/n)
                                  output.n <- n
                                  output <- c(output.lnull, output.n)
                                  names(output) <- c("lnull","n")
                                  return(output)
}

### Calculate Lnull
Lnull.MN <- Lnull.calc.MN(EV1, EV2, EV3, EV4)
Lnull.MN

### Calculate max(R2_CS)
maxR2.MN <- 1 - exp(2*Lnull.MN["lnull"]/Lnull.MN["n"])

### Calculate R2_CS corresponding to an R2_NAGELKERKE of 0.15
R2_CS_NAGEL.MN <- 0.15*maxR2.MN

### Now calculate the N_MN, based on this assumed R2_CS value
### Let sum.p be the total number of predictors considered across all the submodels, which is 16x3
sum.p <- p*3

### Calculate N_MN_NAGEL
N_MN_NAGEL <- sum.p/((S - 1)*log(1 - R2_CS_NAGEL.MN/S))

N_MN_NAGEL
N_DL_NAGEL.2
N_DL_NAGEL.3
N_DL_NAGEL.4



