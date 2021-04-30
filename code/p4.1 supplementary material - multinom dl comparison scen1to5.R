### This program will generate data and fit a multinomial model, and a sequential model
### It will then assess the calibration of this model using both the multinomial recalibration framework,
### and the sequential logistic calibration framework

### It will do this in a simulation, with 1000 simulations per scenario

#install.packages("VGAM")

library(VGAM)
library(foreach)
library(doParallel)

### Set seed
set.seed(1001)

### First define a function that will cerate a dataset according to various beta values
create.dataset <- function(npatt, beta02, beta12, beta22, beta32, beta42, beta52, 
                           beta03, beta13, beta23, beta33, beta43, beta53){
  x1 <- rnorm(npatt, 0, 1)
  x2 <- rnorm(npatt, 0, 1)
  x3 <- rnorm(npatt, 0, 1)
  x4 <- rnorm(npatt, 0, 1)
  x5 <- rnorm(npatt, 0, 1)
  p1 <- 1/(1 + exp(beta02 + x1*beta12 + x2*beta22 + x3*beta32 + x4*beta42 + x5*beta52) + 
             exp(beta03 + x1*beta13 + x2*beta23 + x3*beta33 + x4*beta43 + x5*beta53))
  p2 <- exp(beta02 + x1*beta12 + x2*beta22 + x3*beta32 + x4*beta42 + x5*beta52)/(1 + exp(beta02 + x1*beta12 + x2*beta22 + x3*beta32 + x4*beta42 + x5*beta52) + 
                                                                                   exp(beta03 + x1*beta13 + x2*beta23 + x3*beta33 + x4*beta43 + x5*beta53))
  p3 <- exp(beta03 + x1*beta13 + x2*beta23 + x3*beta33 + x4*beta43 + x5*beta53)/(1 + exp(beta02 + x1*beta12 + x2*beta22 + x3*beta32 + x4*beta42 + x5*beta52) + 
                                                                                   exp(beta03 + x1*beta13 + x2*beta23 + x3*beta33 + x4*beta43 + x5*beta53))
  multinom.probs <- data.frame("p1" = p1, "p2" = p2, "p3" = p3)
  
  ### Now to generate outcome data from a multinomial distribution for each patient
  multinom.outcomes <- matrix(nrow = 3, ncol = npatt)
  for (i in 1:npatt){
    multinom.outcomes[, i] <- rmultinom(n = 1, size = 1, prob = multinom.probs[i, ])
  }
  
  ### Now for each column, need to convert this into a number depending on which row is equal to 1, and add it to a vector
  y.vec <- rep(0, npatt)
  for (i in 1:npatt){
    y.vec[i] <- which(multinom.outcomes[,i] > 0)
  }
  
  
  ### Now combine the outcome data, and the predictor data into a dataset
  devel.data <- data.frame("x1" = x1, "x2" = x2, "x3" = x3, "x4" = x4, "x5" = x5, "y.num" = y.vec)
  
  ### Create a non-numeric vector for y, so it is easier to deal with
  devel.data$y.cat <- rep(0,npatt)
  devel.data$y.cat[devel.data$y.num == 1] <- "cat1"
  devel.data$y.cat[devel.data$y.num == 2] <- "cat2"
  devel.data$y.cat[devel.data$y.num == 3] <- "cat3"
  return(devel.data)
}


### Create ,validation datases
data.valid.s1 <- create.dataset(250000, 0, .5, -0.25, -0.125, 0.25, 0.375,
                                0, 0.375, -0.5, -0.25, -0.375, 0.125)
data.valid.s2 <- create.dataset(250000, -0, .5, -0.25, -0.125, 0.25, 0.375,
                                -0.75, 0.375, -0.5, -0.25, -0.375, 0.125)
data.valid.s3 <- create.dataset(250000, -0.35, .5, -0.25, -0.125, 0.25, 0.375,
                                -0.85, 0.375, -0.5, -0.25, -0.375, 0.125)
data.valid.s4 <- create.dataset(250000, -0.4, .5, -0.25, -0.125, 0.25, 0.375,
                                -1.7, 0.375, -0.5, -0.25, -0.375, 0.125)
data.valid.s5 <- create.dataset(250000, -0.53, .5, -0.25, -0.125, 0.25, 0.375,
                                -2.4, 0.375, -0.5, -0.25, -0.375, 0.125)



### Add the variables that are required for sequential logistic calibration
data.valid.s1$Y1 <- as.numeric(data.valid.s1$y.cat %in% c("cat1"))
data.valid.s1$Y2 <- as.numeric(data.valid.s1$y.cat %in% c("cat2"))
data.valid.s1$Y2[data.valid.s1$Y1 == 1] <- NA
data.valid.s1$Y1 <- factor(data.valid.s1$Y1)
data.valid.s1$Y2 <- factor(data.valid.s1$Y2)

data.valid.s2$Y1 <- as.numeric(data.valid.s2$y.cat %in% c("cat1"))
data.valid.s2$Y2 <- as.numeric(data.valid.s2$y.cat %in% c("cat2"))
data.valid.s2$Y2[data.valid.s2$Y1 == 1] <- NA
data.valid.s2$Y1 <- factor(data.valid.s2$Y1)
data.valid.s2$Y2 <- factor(data.valid.s2$Y2)

data.valid.s3$Y1 <- as.numeric(data.valid.s3$y.cat %in% c("cat1"))
data.valid.s3$Y2 <- as.numeric(data.valid.s3$y.cat %in% c("cat2"))
data.valid.s3$Y2[data.valid.s3$Y1 == 1] <- NA
data.valid.s3$Y1 <- factor(data.valid.s3$Y1)
data.valid.s3$Y2 <- factor(data.valid.s3$Y2)

data.valid.s4$Y1 <- as.numeric(data.valid.s4$y.cat %in% c("cat1"))
data.valid.s4$Y2 <- as.numeric(data.valid.s4$y.cat %in% c("cat2"))
data.valid.s4$Y2[data.valid.s4$Y1 == 1] <- NA
data.valid.s4$Y1 <- factor(data.valid.s4$Y1)
data.valid.s4$Y2 <- factor(data.valid.s4$Y2)

data.valid.s5$Y1 <- as.numeric(data.valid.s5$y.cat %in% c("cat1"))
data.valid.s5$Y2 <- as.numeric(data.valid.s5$y.cat %in% c("cat2"))
data.valid.s5$Y2[data.valid.s5$Y1 == 1] <- NA
data.valid.s5$Y1 <- factor(data.valid.s5$Y1)
data.valid.s5$Y2 <- factor(data.valid.s5$Y2)



### Create datasets in formats that can be used for generation of the linear predictors
lp.data.valid.s1 <- cbind(rep(1,nrow(data.valid.s1)), data.valid.s1$x1, data.valid.s1$x2, data.valid.s1$x3, 
                          data.valid.s1$x4, data.valid.s1$x5)

lp.data.valid.s2 <- cbind(rep(1,nrow(data.valid.s2)), data.valid.s2$x1, data.valid.s2$x2, data.valid.s2$x3, 
                          data.valid.s2$x4, data.valid.s2$x5)

lp.data.valid.s3 <- cbind(rep(1,nrow(data.valid.s3)), data.valid.s3$x1, data.valid.s3$x2, data.valid.s3$x3, 
                          data.valid.s3$x4, data.valid.s3$x5)

lp.data.valid.s4 <- cbind(rep(1,nrow(data.valid.s4)), data.valid.s4$x1, data.valid.s4$x2, data.valid.s4$x3, 
                          data.valid.s4$x4, data.valid.s4$x5)

lp.data.valid.s5 <- cbind(rep(1,nrow(data.valid.s5)), data.valid.s5$x1, data.valid.s5$x2, data.valid.s5$x3, 
                          data.valid.s5$x4, data.valid.s5$x5)




## The above datasets have been set up to allow validation of both multinomial and sequential logistic (NA values have been
## assigned for when Y = 1 for the second sequential logistic model). I also need to add the right outcomes for the distinct logistic 
## regressions, which are both fit to subsets of the validation dataset. Rather than adding multiple outcomes to the same dataset
## and having NA values where the outcomes shouldn't be used (which may get messy), I am just going to create the variables and then
## create seperate datasets where I remove the observations that shouldn't be used
## TL;DR: Creating extra datasets for validation of distinct logistic models

### Scenario 1
## Create the seperate datasets
data.valid.DL1.s1 <- data.valid.s1[(data.valid.s1$y.cat != "cat3"), ]
data.valid.DL2.s1 <- data.valid.s1[(data.valid.s1$y.cat != "cat2"), ]

## Create two outcome variables, the first = 1 if Y = 2, the seonc = 1 if Y = 3
data.valid.DL1.s1$Y.DL1 <- as.numeric(data.valid.DL1.s1$y.cat %in% c("cat2"))
data.valid.DL2.s1$Y.DL2 <- as.numeric(data.valid.DL2.s1$y.cat %in% c("cat3"))

## Create datasets in formats that can be used for generation of the linear predictors
lp.data.valid.DL1.s1 <- cbind(rep(1,nrow(data.valid.DL1.s1)), data.valid.DL1.s1$x1, data.valid.DL1.s1$x2, data.valid.DL1.s1$x3, 
                              data.valid.DL1.s1$x4, data.valid.DL1.s1$x5)
lp.data.valid.DL2.s1 <- cbind(rep(1,nrow(data.valid.DL2.s1)), data.valid.DL2.s1$x1, data.valid.DL2.s1$x2, data.valid.DL2.s1$x3, 
                              data.valid.DL2.s1$x4, data.valid.DL2.s1$x5)


### Scenario 2
## Create the seperate datasets
data.valid.DL1.s2 <- data.valid.s2[(data.valid.s2$y.cat != "cat3"), ]
data.valid.DL2.s2 <- data.valid.s2[(data.valid.s2$y.cat != "cat2"), ]

## Create two outcome variables, the first = 1 if Y = 2, the seonc = 1 if Y = 3
data.valid.DL1.s2$Y.DL1 <- as.numeric(data.valid.DL1.s2$y.cat %in% c("cat2"))
data.valid.DL2.s2$Y.DL2 <- as.numeric(data.valid.DL2.s2$y.cat %in% c("cat3"))

## Create datasets in formats that can be used for generation of the linear predictors
lp.data.valid.DL1.s2 <- cbind(rep(1,nrow(data.valid.DL1.s2)), data.valid.DL1.s2$x1, data.valid.DL1.s2$x2, data.valid.DL1.s2$x3, 
                              data.valid.DL1.s2$x4, data.valid.DL1.s2$x5)
lp.data.valid.DL2.s2 <- cbind(rep(1,nrow(data.valid.DL2.s2)), data.valid.DL2.s2$x1, data.valid.DL2.s2$x2, data.valid.DL2.s2$x3, 
                              data.valid.DL2.s2$x4, data.valid.DL2.s2$x5)


### Scenario 3
## Create the seperate datasets
data.valid.DL1.s3 <- data.valid.s3[(data.valid.s3$y.cat != "cat3"), ]
data.valid.DL2.s3 <- data.valid.s3[(data.valid.s3$y.cat != "cat2"), ]

## Create two outcome variables, the first = 1 if Y = 2, the seonc = 1 if Y = 3
data.valid.DL1.s3$Y.DL1 <- as.numeric(data.valid.DL1.s3$y.cat %in% c("cat2"))
data.valid.DL2.s3$Y.DL2 <- as.numeric(data.valid.DL2.s3$y.cat %in% c("cat3"))

## Create datasets in formats that can be used for generation of the linear predictors
lp.data.valid.DL1.s3 <- cbind(rep(1,nrow(data.valid.DL1.s3)), data.valid.DL1.s3$x1, data.valid.DL1.s3$x2, data.valid.DL1.s3$x3, 
                              data.valid.DL1.s3$x4, data.valid.DL1.s3$x5)
lp.data.valid.DL2.s3 <- cbind(rep(1,nrow(data.valid.DL2.s3)), data.valid.DL2.s3$x1, data.valid.DL2.s3$x2, data.valid.DL2.s3$x3, 
                              data.valid.DL2.s3$x4, data.valid.DL2.s3$x5)


### Scenario 4
## Create the seperate datasets
data.valid.DL1.s4 <- data.valid.s4[(data.valid.s4$y.cat != "cat3"), ]
data.valid.DL2.s4 <- data.valid.s4[(data.valid.s4$y.cat != "cat2"), ]

## Create two outcome variables, the first = 1 if Y = 2, the seonc = 1 if Y = 3
data.valid.DL1.s4$Y.DL1 <- as.numeric(data.valid.DL1.s4$y.cat %in% c("cat2"))
data.valid.DL2.s4$Y.DL2 <- as.numeric(data.valid.DL2.s4$y.cat %in% c("cat3"))

## Create datasets in formats that can be used for generation of the linear predictors
lp.data.valid.DL1.s4 <- cbind(rep(1,nrow(data.valid.DL1.s4)), data.valid.DL1.s4$x1, data.valid.DL1.s4$x2, data.valid.DL1.s4$x3, 
                              data.valid.DL1.s4$x4, data.valid.DL1.s4$x5)
lp.data.valid.DL2.s4 <- cbind(rep(1,nrow(data.valid.DL2.s4)), data.valid.DL2.s4$x1, data.valid.DL2.s4$x2, data.valid.DL2.s4$x3, 
                              data.valid.DL2.s4$x4, data.valid.DL2.s4$x5)


### Scenario 5
## Create the seperate datasets
data.valid.DL1.s5 <- data.valid.s5[(data.valid.s5$y.cat != "cat3"), ]
data.valid.DL2.s5 <- data.valid.s5[(data.valid.s5$y.cat != "cat2"), ]

## Create two outcome variables, the first = 1 if Y = 2, the seonc = 1 if Y = 3
data.valid.DL1.s5$Y.DL1 <- as.numeric(data.valid.DL1.s5$y.cat %in% c("cat2"))
data.valid.DL2.s5$Y.DL2 <- as.numeric(data.valid.DL2.s5$y.cat %in% c("cat3"))

## Create datasets in formats that can be used for generation of the linear predictors
lp.data.valid.DL1.s5 <- cbind(rep(1,nrow(data.valid.DL1.s5)), data.valid.DL1.s5$x1, data.valid.DL1.s5$x2, data.valid.DL1.s5$x3, 
                              data.valid.DL1.s5$x4, data.valid.DL1.s5$x5)
lp.data.valid.DL2.s5 <- cbind(rep(1,nrow(data.valid.DL2.s5)), data.valid.DL2.s5$x1, data.valid.DL2.s5$x2, data.valid.DL2.s5$x3, 
                              data.valid.DL2.s5$x4, data.valid.DL2.s5$x5)


#####################################################################################################################
### Write a function to fit multinomial, sequential logistic and distinct logistic models in the development dataset
### and then test their calibration in the validation dataset, and also calculate heuristic shrinkage factors of each model
#####################################################################################################################
get.slopes.loglik <- function(data.devel.in, data.valid.in, lp.data.valid.in, data.valid.DL1.in, lp.data.valid.DL1.in, 
                              data.valid.DL2.in, lp.data.valid.DL2.in){
  

  ##########################################
  ### Develop a multinoial logistic model
  ##########################################
  
  ### Want to develop a model in devel
  multinom.model <- vglm(y.cat ~ x1 + x2 + x3 + x4 + x5, family = multinomial(ref = "cat1"), data = data.devel.in)
  
  ## Extract coefficients
  coeffs.multinom.1 <- multinom.model@coefficients[c(1,3,5,7,9,11)]
  coeffs.multinom.2 <- multinom.model@coefficients[c(2,4,6,8,10,12)]
  
  ## Calculate linear predictors
  lp.valid.multinom.1 <- lp.data.valid.in %*% coeffs.multinom.1
  lp.valid.multinom.2 <- lp.data.valid.in %*% coeffs.multinom.2
  
  ## Also store the standard errors of the effect estimates
  se.multinom.lp1 <- coef(summary(multinom.model))[c(1,3,5,7,9,11), c("Std. Error")]
  se.multinom.lp2 <- coef(summary(multinom.model))[c(2,4,6,8,10,12), c("Std. Error")]
  
  
  ##########################################
  ### Develop distinct logistic models
  ##########################################
  
  ## Now need to create two other datasets, retaining only the observations pertinant to the two distinct logistic models
  data.devel.in.DL1 <- data.devel.in[(data.devel.in$y.cat != "cat3"), ]
  data.devel.in.DL2 <- data.devel.in[(data.devel.in$y.cat != "cat2"), ]
  
  ## First create two outcome variables, the first = 1 if Y = 2, the seonc = 1 if Y = 3
  data.devel.in.DL1$Y.DL1 <- as.numeric(data.devel.in.DL1$y.cat %in% c("cat2"))
  data.devel.in.DL2$Y.DL2 <- as.numeric(data.devel.in.DL2$y.cat %in% c("cat3"))
 
  
  ## Now need to fit both the models
  dislog.model.1 <- glm(Y.DL1 ~ x1 + x2 + x3 + x4 + x5, family = binomial(link = "logit"), 
                  data = data.devel.in.DL1)
  dislog.model.2 <- glm(Y.DL2 ~ x1 + x2 + x3 + x4 + x5, family = binomial(link = "logit"), 
                  data = data.devel.in.DL2)
  
  ## Extract coefficients
  coeffs.dislog.1 <- dislog.model.1$coefficients
  coeffs.dislog.2 <- dislog.model.2$coefficients
  
  ## Also extract standard errors
  se.dl.lp1 <- coef(summary(dislog.model.1))[, c("Std. Error")]
  se.dl.lp2 <- coef(summary(dislog.model.2))[, c("Std. Error")]
  
  ## Calculate linear predictors in validaiton dataset
  lp.valid.dislog.1 <- lp.data.valid.DL1.in %*% coeffs.dislog.1
  lp.valid.dislog.2 <- lp.data.valid.DL2.in %*% coeffs.dislog.2
  
  ## Also calculate linear predictor's for the entire valdation dataset, to be able to
  ## do the multinomial validation of the DL's
  lp.valid.all.dislog.1 <- lp.data.valid.in %*% coeffs.dislog.1
  lp.valid.all.dislog.2 <- lp.data.valid.in %*% coeffs.dislog.2
  
  
  
  ########################################################
  ### Now to calculate the calibration slopes
  ### This will be done using either a multinomial recalibration framework, or just standard
  ### calibration intercept/slope for the sequential and distinct logistic regressions
  ### Note for the multinomial and sequential logistic, I also have written code to calibrate the risk scores
  ### using the alternative calibration method (i.e. generate risks using seq log, then calibrate using multinomial, and vice versa)
  ### This was done experimentally and the code is commented out
  ########################################################
  
  
  ########################################
  ### Model developed: Multinomial 
  ### Calibration: Multinomial
  ########################################
  
  ## Create a dataset for the recalibration model to be fit on
  recal.data <- data.frame("y.cat" = data.valid.in$y.cat,
                           "lp.1" = lp.valid.multinom.1, 
                           "lp.2" = lp.valid.multinom.2)
  
  ## Now do the recalibration
  i <- diag(2)
  i2 <- rbind(1, 0)
  i3 <- rbind(0, 1)
  clist <- list("(Intercept)" = i, "lp.1" = i2, "lp.2" = i3)
  clist
  
  ## Fit the recalibration model (contains intercepts and slopes)
  devel.m.cali.m <- vgam(y.cat ~ lp.1 + lp.2, family = multinomial(ref = "cat1"), constraints = clist, data = recal.data)
  
  
  
  ########################################
  ### Model developed: Distinct logistic 
  ### Calibration: Multinomial
  ########################################
  
  ## Create a dataset for the recalibration model to be fit on
  recal.data <- data.frame("y.cat" = data.valid.in$y.cat,
                           "lp.1" = lp.valid.all.dislog.1, 
                           "lp.2" = lp.valid.all.dislog.2)
  
  ## Now do the recalibration
  i <- diag(2)
  i2 <- rbind(1, 0)
  i3 <- rbind(0, 1)
  clist <- list("(Intercept)" = i, "lp.1" = i2, "lp.2" = i3)
  clist
  
  ## Fit the recalibration model (contains intercepts and slopes)
  devel.d.cali.m <- vgam(y.cat ~ lp.1 + lp.2, family = multinomial(ref = "cat1"), constraints = clist, data = recal.data)
  
  
  ### Create output object
  output.object <- vector("list", 4)
  names(output.object) <- c("devel.m.cali.m", "devel.d.cali.m","se.multinom","se.dl")
  
  output.object$devel.m.cali.m <- devel.m.cali.m@coefficients
  
  output.object$devel.d.cali.m <- devel.d.cali.m@coefficients
  
  output.object$se.multinom <- c(se.multinom.lp1, se.multinom.lp2)
    
  output.object$se.dl <- c(se.dl.lp1, se.dl.lp2)
  
  return(output.object)}



### Now for each scenario, I want to run a simulation for a variety of sample size
### I will write a function that will do this

## Define the function
## Define the function
run.multinom.sim <- function(a.beta02, a.beta12, a.beta22, a.beta32, a.beta42, a.beta52, 
                             a.beta03, a.beta13, a.beta23, a.beta33, a.beta43, a.beta53, 
                             data.valid.inn, lp.data.valid.inn, data.valid.DL1.inn, lp.data.valid.DL1.inn, 
                             data.valid.DL2.inn, lp.data.valid.DL2.inn, n.sim, sample.size){
  # Important to use the correct validation data, corresponding to the beta values that were inputted
  # a. prefix is just to distinguish from betaX defined in previous functions, to help with any debugging
     
    ## Create an output dataset for the sample size
    output.data <- matrix(nrow = n.sim, ncol = 32)
    
    ## Give the columns appropiate names
    colnames(output.data) <- c("Int.devel.m.cali.m.1","Int.devel.m.cali.m.2","Slope.devel.m.cali.m.1","Slope.devel.m.cali.m.2",
                               "Int.devel.d.cali.m.1","Int.devel.d.cali.m.2","Slope.devel.d.cali.m.1","Slope.devel.d.cali.m.2",
                               "se.m.1.1","se.m.1.2","se.m.1.3","se.m.1.4","se.m.1.5","se.m.1.6",
                               "se.m.2.1","se.m.2.2","se.m.2.3","se.m.2.4","se.m.2.5","se.m.2.6",
                               "se.dl.1.1","se.dl.1.2","se.dl.1.3","se.dl.1.4","se.dl.1.5","se.dl.1.6",
                               "se.dl.2.1","se.dl.2.2","se.dl.2.3","se.dl.2.4","se.dl.2.5","se.dl.2.6")
    
    ## Run through the simulation, 1000 times
    for (j in 1:n.sim){
      
      ## Create dataset
      data.devell <- create.dataset(sample.size, a.beta02, a.beta12, a.beta22, a.beta32, a.beta42, a.beta52, 
                                    a.beta03, a.beta13, a.beta23, a.beta33, a.beta43, a.beta53)
      
      ## Calculate the calibration slopes and intercepts in the validation dataset, and the loglikelihood of the
      ## model
      slopes.loglik <- get.slopes.loglik(data.devell, data.valid.inn, lp.data.valid.inn, data.valid.DL1.inn, lp.data.valid.DL1.inn, 
                                         data.valid.DL2.inn, lp.data.valid.DL2.inn)
      output.data[j, 1:4] <- slopes.loglik$devel.m.cali.m
      output.data[j, 5:8] <- slopes.loglik$devel.d.cali.m
      output.data[j, 9:20] <- slopes.loglik$se.multinom
      output.data[j, 21:32] <- slopes.loglik$se.dl}
    
    print(warnings())
    return(output.data)
  }


### Scenario 1
print("start")
cl <- makeCluster(9)
registerDoParallel(9)
scenario1.res<-(foreach(input=c(100,250,500,1000,540,564,576,572,606), .combine=list, .multicombine=TRUE)
                      %dopar%{run.multinom.sim(0, .5, -0.25, -0.125, 0.25, 0.375,
                                               0, 0.375, -0.5, -0.25, -0.375, 0.125, 
                                               data.valid.s1, lp.data.valid.s1, data.valid.DL1.s1, lp.data.valid.DL1.s1, 
                                               data.valid.DL2.s1, lp.data.valid.DL2.s1, 1000, input)
                      })
stopCluster(cl)

print("scenario 1 done")
#warnings()
Sys.time()
save.image("R_out/multinom dl comparison scen1to5.RData")


## Scenario 2
## Define sample sizes
cl <- makeCluster(9)
registerDoParallel(9)
scenario2.res<-(foreach(input=c(100,250,500,1000,569,699,628,570,571), .combine=list, .multicombine=TRUE)
                %dopar%{run.multinom.sim(-0, .5, -0.25, -0.125, 0.25, 0.375,
                                         -0.75, 0.375, -0.5, -0.25, -0.375, 0.125, 
                                         data.valid.s2, lp.data.valid.s2, data.valid.DL1.s2, lp.data.valid.DL1.s2, 
                                         data.valid.DL2.s2, lp.data.valid.DL2.s2, 1000,input)
                })
stopCluster(cl)



print("scenario 2 done")
#warnings()
Sys.time()
save.image("R_out/multinom dl comparison scen1to5.RData")


### Scenario 3
## Define sample sizes
cl <- makeCluster(9)
registerDoParallel(9)
scenario3.res<-(foreach(input=c(100,250,500,1000,565,714,582,565,566), .combine=list, .multicombine=TRUE)
                %dopar%{run.multinom.sim(-0.35, .5, -0.25, -0.125, 0.25, 0.375,
                                         -0.85, 0.375, -0.5, -0.25, -0.375, 0.125, 
                                         data.valid.s3, lp.data.valid.s3, data.valid.DL1.s3, lp.data.valid.DL1.s3, 
                                         data.valid.DL2.s3, lp.data.valid.DL2.s3, 1000,input)
                })
stopCluster(cl)


print("scenario 3 done")
#warnings()
Sys.time()
save.image("R_out/multinom dl comparison scen1to5.RData")


### Scenario 4
## Define sample sizes
cl <- makeCluster(9)
registerDoParallel(9)
scenario4.res<-(foreach(input=c(250,500,1000,647,1105,901,813,1025), .combine=list, .multicombine=TRUE)
                %dopar%{run.multinom.sim(-0.4, .5, -0.25, -0.125, 0.25, 0.375,
                                         -1.7, 0.375, -0.5, -0.25, -0.375, 0.125, 
                                         data.valid.s4, lp.data.valid.s4, data.valid.DL1.s4, lp.data.valid.DL1.s4, 
                                         data.valid.DL2.s4, lp.data.valid.DL2.s4, 1000,input)
                })
stopCluster(cl)


print("scenario 4 done")
#warnings()
Sys.time()
save.image("R_out/multinom dl comparison scen1to5.RData")



### Scenario 5
## Define sample sizes
cl <- makeCluster(9)
registerDoParallel(9)
scenario5.res<-(foreach(input=c(250,500,1000,706,1742,1457,1156,2003), .combine=list, .multicombine=TRUE)
                %dopar%{run.multinom.sim(-0.53, .5, -0.25, -0.125, 0.25, 0.375,
                                         -2.4, 0.375, -0.5, -0.25, -0.375, 0.125, 
                                         data.valid.s5, lp.data.valid.s5, data.valid.DL1.s5, lp.data.valid.DL1.s5, 
                                         data.valid.DL2.s5, lp.data.valid.DL2.s5, 1000,input)
                })
stopCluster(cl)


print("scenario 5 done")
#warnings()
Sys.time()
save.image("R_out/multinom dl comparison scen1to5.RData")

