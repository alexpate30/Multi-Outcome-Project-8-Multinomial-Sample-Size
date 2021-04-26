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


### Create validation datasets
data.valid.s1 <- create.dataset(500000, 0, .5, -0.25, -0.125, 0.25, 0.375,
                                0, 0.375, -0.5, -0.25, -0.375, 0.125)
data.valid.s2 <- create.dataset(500000, -0, .5, -0.25, -0.125, 0.25, 0.375,
                                -0.75, 0.375, -0.5, -0.25, -0.375, 0.125)
data.valid.s3 <- create.dataset(500000, -0.35, .5, -0.25, -0.125, 0.25, 0.375,
                                -0.85, 0.375, -0.5, -0.25, -0.375, 0.125)
data.valid.s4 <- create.dataset(500000, -0.4, .5, -0.25, -0.125, 0.25, 0.375,
                                -1.7, 0.375, -0.5, -0.25, -0.375, 0.125)
data.valid.s5 <- create.dataset(500000, -0.53, .5, -0.25, -0.125, 0.25, 0.375,
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
  multinom.model <- vgam(y.cat ~ x1 + x2 + x3 + x4 + x5, family = multinomial(ref = "cat1"), data = data.devel.in)
  
  ## Extract coefficients
  coeffs.multinom.1 <- multinom.model@coefficients[c(1,3,5,7,9,11)]
  coeffs.multinom.2 <- multinom.model@coefficients[c(2,4,6,8,10,12)]
  
  ## Calculate linear predictors
  lp.valid.multinom.1 <- lp.data.valid.in %*% coeffs.multinom.1
  lp.valid.multinom.2 <- lp.data.valid.in %*% coeffs.multinom.2
  
  
  ##########################################
  ### Develop a sequential logistic model
  ##########################################
  
  ### Create dummy variables in the development dataset
  data.devel.in$Y1 <- as.numeric(data.devel.in$y.cat %in% c("cat1"))
  data.devel.in$Y2 <- as.numeric(data.devel.in$y.cat %in% c("cat2"))
  data.devel.in$Y2[data.devel.in$Y1 == 1] <- NA
  
  data.devel.in$Y1 <- factor(data.devel.in$Y1)
  data.devel.in$Y2 <- factor(data.devel.in$Y2)
  
  ### Fit the models
  seqlog.model.1 <- glm(Y1 ~ x1 + x2 + x3 + x4 + x5, family = binomial(link = "logit"), 
                        data = data.devel.in)
  seqlog.model.2 <- glm(Y2 ~ x1 + x2 + x3 + x4 + x5, family = binomial(link = "logit"), 
                        data = data.devel.in[!is.na(data.devel.in$Y2), ])
  
  ## Extract coefficients
  coeffs.seqlog.1 <- seqlog.model.1$coefficients
  coeffs.seqlog.2 <- seqlog.model.2$coefficients
  
  ## Calculate linear predictors
  lp.valid.seqlog.1 <- lp.data.valid.in %*% coeffs.seqlog.1
  lp.valid.seqlog.2 <- lp.data.valid.in %*% coeffs.seqlog.2
  
  
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
  
  ## Calculate linear predictors in validaiton dataset
  lp.valid.dislog.1 <- lp.data.valid.DL1.in %*% coeffs.dislog.1
  lp.valid.dislog.2 <- lp.data.valid.DL2.in %*% coeffs.dislog.2
  
  
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
  ### Model developed: Sequential logistic 
  ### Calibration: Sequential logistic
  ########################################
  
  ## Create a dataset for the recalibration model to be fit on
  recal.data <- data.frame("Y1" = data.valid.in$Y1,
                           "Y2" = data.valid.in$Y2,
                           "lp.1" = lp.valid.seqlog.1, 
                           "lp.2" = lp.valid.seqlog.2)
  
  
  ## Do the calibration
  devel.s.cali.s.1 <- glm(Y1 ~ lp.1, family = binomial(link = "logit"), 
                          data = recal.data)
  devel.s.cali.s.2 <- glm(Y2 ~ lp.2, family = binomial(link = "logit"), 
                          data = recal.data[!is.na(recal.data$Y2), ])
  
  
  ########################################
  ### Model developed: Distinct logistic 
  ### Calibration: Distinct logistic
  ########################################

  ## Create datasets for the recalibration model to be fit on
  recal.data.DL1 <- data.frame("Y.DL1" = data.valid.DL1.in$Y.DL1,
                           "lp.1" = lp.valid.dislog.1)
  recal.data.DL2 <- data.frame("Y.DL2" = data.valid.DL2.in$Y.DL2,
                               "lp.2" = lp.valid.dislog.2)
  
  ## Do the calibration
  devel.d.cali.d.1 <- glm(Y.DL1 ~ lp.1, family = binomial(link = "logit"), 
                          data = recal.data.DL1)
  devel.d.cali.d.2 <- glm(Y.DL2 ~ lp.2, family = binomial(link = "logit"), 
                          data = recal.data.DL2)
  
  
  ########################################
  ### Model developed: Multinomial
  ### Calibration: Sequential logistic
  ########################################
  
#   ## Create a dataset for the recalibration model to be fit on
#   ## Given we are using the sequential logistic recalibration, we need those outcomes
#   ## However the predictors must be from the multinomial
#   recal.data <- data.frame("Y1" = data.valid.in$Y1,
#                            "Y2" = data.valid.in$Y2,
#                            "lp.1" = lp.valid.multinom.1, 
#                            "lp.2" = lp.valid.multinom.2)
#   
#   ## We must then generate the required entities for the sequential logistic calibration
#   ## These are log(P(Y=1)/(1-P(Y=1))), and log(P(Y=2)/(P(Y=3)))
#   ## These are the two regressions that are actually done in the sequential logistic
#   ## These can be calculated from the linear predictors from the multinomial model
#   
#   ##First calculate the probabilities of having each outcome
#   recal.data$risk.y1 <- 1/(1 + exp(recal.data$lp.1) + exp(recal.data$lp.2))
#   recal.data$risk.y2 <- exp(recal.data$lp.1)/(1 + exp(recal.data$lp.1) + exp(recal.data$lp.2))
#   recal.data$risk.y3 <- exp(recal.data$lp.2)/(1 + exp(recal.data$lp.1) + exp(recal.data$lp.2))
#   
#   ## Now derive what will be the predictor for the two calibration models
#   recal.data$mod1pred <- log(recal.data$risk.y1/(1-recal.data$risk.y1))
#   recal.data$mod2pred <- log(recal.data$risk.y2/recal.data$risk.y3)
#   
#   ## Now do the calibration
#   devel.m.cali.s.1 <- glm(Y1 ~ mod1pred, family = binomial(link = "logit"), 
#                           data = recal.data)
#   devel.m.cali.s.2 <- glm(Y2 ~ mod2pred, family = binomial(link = "logit"), 
#                           data = recal.data[!is.na(recal.data$Y2), ])
#   
  
  
  ########################################
  ### Model developed: Sequential logistic
  ### Calibration: Multinomial
  ########################################
  
#   ## Create a dataset for the recalibration model to be fit on
#   ## Given we are using the sequential logisticmultinomial recalibration, we need those outcomes
#   ## However the predictors must be from the sequential logistic
#   recal.data <- data.frame("y.cat" = data.valid.in$y.cat,
#                            "lp.1" = lp.valid.seqlog.1, 
#                            "lp.2" = lp.valid.seqlog.2)
#   
#   
#   ## First calculate probabilities of each outcome
#   recal.data$risk.y1 <- exp(recal.data$lp.1)/(1 + exp(recal.data$lp.1))
#   recal.data$risk.y2 <- (1/(1+exp(-recal.data$lp.2)))*(1-recal.data$risk.y1)
#   recal.data$risk.y3 <- (1/(1+exp(recal.data$lp.2)))*(1-recal.data$risk.y1)
#   
#   
#   ## Now derive what will be the predictor for the two calibration models
#   recal.data$mod1pred <- log(recal.data$risk.y2/recal.data$risk.y1)
#   recal.data$mod2pred <- log(recal.data$risk.y3/recal.data$risk.y1)
#   
#   
#   ## Now do the recalibration
#   i <- diag(2)
#   i2 <- rbind(1, 0)
#   i3 <- rbind(0, 1)
#   clist <- list("(Intercept)" = i, "mod1pred" = i2, "mod2pred" = i3)
#   clist
#   
#   ## Fit the recalibration model (contains intercepts and slopes)
#   devel.s.cali.m <- vgam(y.cat ~ mod1pred + mod2pred, family = multinomial(ref = "cat1"), 
#                          constraints = clist, data = recal.data)
#   
  
  
  ### I have now developed models using both approaches, and tested calibration in the validatiobn
  ### dataset using both approaches. Finally want to record the heuristic shrinkage factors of each 
  ### model. For the multinomial model, we just have one 'overall' heuristic shrinkage factor. For
  ### the sequential logistic, there are two shrinkage factors (one for each sequential logistic model).
  
  ### Note that the heuristic shrinkage factor applies to the fitted model only, and therefore this process
  ### will not need to be repeated based on which recalibration framework is used.
  
  ### Finally, I will actually just output the likelihood ratios, as then I can make different decisions about
  ### what p should be when summarising the data.
  
  
  #########################################
  ### Likelihood ratio of multinomial model
  #########################################
  
  ### First need to fit an intercept only model to the development dataset
  multinom.model.null <- vgam(y.cat ~ 1, family = multinomial(ref = "cat1"), data = data.devel.in)
  
  ### Then calculate the likelihood ratio
  multinom.LR <- -2*(multinom.model.null@criterion$loglikelihood - multinom.model@criterion$loglikelihood) 
  
  
  ###############################################################################################
  ### Likelihood ratio of each categories contribution to the likelihood of the multinomial model
  ###############################################################################################
  
  #### Want to calculate the likelihood ratios of the contributions by the specific outcomes
  #### To do this, need to calculate the contributions seperately from the full model and null model
  
  ##########################
  #### Start with full model
  ##########################
  
  ## Have to start by calculating the risk scores for each individual in the development cohort
  ## Start by creating a dataset used to generate the linear predictor
  lp.data.devel <- cbind(rep(1,nrow(data.devel.in)), data.devel.in$x1, data.devel.in$x2, data.devel.in$x3, 
                         data.devel.in$x4, data.devel.in$x5)
  
  ## Calculate the linear predictors
  lp.devel.multinom.1 <- lp.data.devel %*% coeffs.multinom.1
  lp.devel.multinom.2 <- lp.data.devel %*% coeffs.multinom.2
  
  ## Now put these into a data frame for ease of use
  data.devel.risks <- data.frame("y.cat" = data.devel.in$y.cat,
                                 "lp.devel.y2" = lp.devel.multinom.1, 
                                 "lp.devel.y3" = lp.devel.multinom.2)
  
  ## Now use these to calculate risk scores
  data.devel.risks$risk.devel.y1 <- 1/(1 + exp(data.devel.risks$lp.devel.y2) + exp(data.devel.risks$lp.devel.y3))
  data.devel.risks$risk.devel.y2 <- exp(data.devel.risks$lp.devel.y2)/(1 + exp(data.devel.risks$lp.devel.y2) + exp(data.devel.risks$lp.devel.y3))
  data.devel.risks$risk.devel.y3 <- exp(data.devel.risks$lp.devel.y3)/(1 + exp(data.devel.risks$lp.devel.y2) + exp(data.devel.risks$lp.devel.y3))
  
  
  ## Now to calculate the loglikihood contributed by each outcome category
  loglik.y1 <- sum(as.numeric(data.devel.risks$y.cat == "cat1")*log(data.devel.risks$risk.devel.y1))
  loglik.y2 <- sum(as.numeric(data.devel.risks$y.cat == "cat2")*log(data.devel.risks$risk.devel.y2))
  loglik.y3 <- sum(as.numeric(data.devel.risks$y.cat == "cat3")*log(data.devel.risks$risk.devel.y3))
  
  
  ################################
  #### Do the same with null model
  ################################
  
  ## The linear predictor for every patient is the same, and is equal to the value of the intercept
  ## Put these into a data frame
  data.devel.risks.null <- data.frame("y.cat" = data.devel.in$y.cat,
                                      "lp.null.devel.y2" = rep(multinom.model.null@coefficients[1], nrow(data.devel.in)), 
                                      "lp.null.devel.y3" = rep(multinom.model.null@coefficients[2], nrow(data.devel.in)))
  
  ## Now use these to calculate risk scores
  data.devel.risks.null$risk.null.devel.y1 <- 1/(1 + exp(data.devel.risks.null$lp.null.devel.y2) + exp(data.devel.risks.null$lp.null.devel.y3))
  data.devel.risks.null$risk.null.devel.y2 <- exp(data.devel.risks.null$lp.null.devel.y2)/(1 + exp(data.devel.risks.null$lp.null.devel.y2) + exp(data.devel.risks.null$lp.null.devel.y3))
  data.devel.risks.null$risk.null.devel.y3 <- exp(data.devel.risks.null$lp.null.devel.y3)/(1 + exp(data.devel.risks.null$lp.null.devel.y2) + exp(data.devel.risks.null$lp.null.devel.y3))
  
  ## Now to calculate the loglikihood
  loglik.null.y1 <- sum(as.numeric(data.devel.risks.null$y.cat == "cat1")*log(data.devel.risks.null$risk.null.devel.y1))
  loglik.null.y2 <- sum(as.numeric(data.devel.risks.null$y.cat == "cat2")*log(data.devel.risks.null$risk.null.devel.y2))
  loglik.null.y3 <- sum(as.numeric(data.devel.risks.null$y.cat == "cat3")*log(data.devel.risks.null$risk.null.devel.y3))
  
  ### Finally calculate the likelihood ratios from these
  multinom.LR.1 <- -2*(loglik.null.y1 - loglik.y1) 
  multinom.LR.2 <- -2*(loglik.null.y2 - loglik.y2) 
  multinom.LR.3 <- -2*(loglik.null.y3 - loglik.y3) 
  
  ## Also save the loglikelihood contributions to the full model, in case these are of interest
  multinom.LL.1 <- loglik.y1
  multinom.LL.2 <- loglik.y2
  multinom.LL.3 <- loglik.y3
  

  ###################################################
  ### Likelihood ratios of sequential logistic models
  ###################################################
  
  ### First need to fit an intercept only model to the development dataset
  seqlog.model.null.1 <- glm(Y1 ~ 1, family = binomial(link = "logit"), 
                             data = data.devel.in)
  seqlog.model.null.2 <- glm(Y2 ~ 1, family = binomial(link = "logit"), 
                             data = data.devel.in[!is.na(data.devel.in$Y2), ])
  
  ### Then calculate the likelihood ratio for each model
  seqlog.LR.1 <- -2*(logLik(seqlog.model.null.1) - logLik(seqlog.model.1))
  seqlog.LR.2 <- -2*(logLik(seqlog.model.null.2) - logLik(seqlog.model.2))
  

  #################################################
  ### Likelihood ratios of distinct logistic models
  #################################################

  ## Fit intercept only models
  dislog.model.null.1 <- glm(Y.DL1 ~ 1, family = binomial(link = "logit"), 
                        data = data.devel.in.DL1)
  dislog.model.null.2 <- glm(Y.DL2 ~ 1, family = binomial(link = "logit"), 
                        data = data.devel.in.DL2)

  ### Then calculate the likelihood ratio for each model
  dislog.LR.1 <- -2*(logLik(dislog.model.null.1) - logLik(dislog.model.1))
  dislog.LR.2 <- -2*(logLik(dislog.model.null.2) - logLik(dislog.model.2))


  #################################################
  ### Put all this into an object for output
  #################################################
  
  #   devel.m.cali.m
  #   devel.s.cali.m
  #   
  #   devel.m.cali.s.1
  #   devel.m.cali.s.2
  #   
  #   devel.s.cali.s.1
  #   devel.s.cali.s.2
  
  ### Create output object
  output.object <- vector("list", 6)
  names(output.object) <- c("devel.m.cali.m", "devel.m.cali.s","devel.s.cali.m","devel.s.cali.s","devel.d.cali.d","LR")
  
  output.object$devel.m.cali.m <- devel.m.cali.m@coefficients
  
#   output.object$devel.m.cali.s <- c(devel.m.cali.s.1$coefficients[1],devel.m.cali.s.2$coefficients[1],
#                                     devel.m.cali.s.1$coefficients[2],devel.m.cali.s.2$coefficients[2])
#   
#   output.object$devel.s.cali.m <- devel.s.cali.m@coefficients
  
  output.object$devel.s.cali.s <- c(devel.s.cali.s.1$coefficients[1],devel.s.cali.s.2$coefficients[1],
                                    devel.s.cali.s.1$coefficients[2],devel.s.cali.s.2$coefficients[2])

  output.object$devel.d.cali.d <- c(devel.d.cali.d.1$coefficients[1],devel.d.cali.d.2$coefficients[1],
                                  devel.d.cali.d.1$coefficients[2],devel.d.cali.d.2$coefficients[2])
  
  output.object$LR <- c(multinom.LR, seqlog.LR.1, seqlog.LR.2, dislog.LR.1, dislog.LR.2, multinom.LR.1, multinom.LR.2, multinom.LR.3, 
                        multinom.LL.1, multinom.LL.2, multinom.LL.3)
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
    output.data <- matrix(nrow = n.sim, ncol = 31)
    
    ## Give the columns appropiate names
    colnames(output.data) <- c("Int.devel.m.cali.m.1","Int.devel.m.cali.m.2","Slope.devel.m.cali.m.1","Slope.devel.m.cali.m.2",
                               "Int.devel.m.cali.s.1","Int.devel.m.cali.s.2","Slope.devel.m.cali.s.1","Slope.devel.m.cali.s.2",
                               "Int.devel.s.cali.m.1","Int.devel.s.cali.m.2","Slope.devel.s.cali.m.1","Slope.devel.s.cali.m.2",
                               "Int.devel.s.cali.s.1","Int.devel.s.cali.s.2","Slope.devel.s.cali.s.1","Slope.devel.s.cali.s.2",
                               "Int.devel.d.cali.d.1","Int.devel.d.cali.d.2","Slope.devel.d.cali.d.1","Slope.devel.d.cali.d.2",
                               "LR.multinom", "LR.seqlog.1", "LR.seqlog.2","LR.dislog.1", "LR.dislog.2", "LR.multinom.1", "LR.multinom.2", "LR.multinom.3",
                               "LL.multinom.1", "LL.multinom.2", "LL.multinom.3")
    
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
      #output.data[j, 5:8] <- slopes.loglik$devel.m.cali.s
      #output.data[j, 9:12] <- slopes.loglik$devel.s.cali.m
      output.data[j, 13:16] <- slopes.loglik$devel.s.cali.s
      output.data[j, 17:20] <- slopes.loglik$devel.d.cali.d
      output.data[j, 21:31] <- slopes.loglik$LR}
    
    print(warnings())
    return(output.data)
  }


### Scenario 1
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
warnings()
Sys.time()
save.image("R_out/multinom simulation scen1to5 500000val.RData")



### Scenario 2
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
warnings()
Sys.time()
save.image("R_out/multinom simulation scen1to5 500000val.RData")


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
warnings()
Sys.time()
save.image("R_out/multinom simulation scen1to5 500000val.RData")


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
warnings()
Sys.time()
save.image("R_out/multinom simulation scen1to5 500000val.RData")



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
warnings()
Sys.time()
save.image("R_out/multinom simulation scen1to5 500000val.RData")

