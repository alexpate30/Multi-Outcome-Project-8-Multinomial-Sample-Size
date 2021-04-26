##### This R code will estimate R2_CS for each distinct logistic regression model, following the process of Riley et al (doi: 10.1002/sim.8806)
##### The C statistic values on which the calculation will be based, are extracted from van Calster et al (doi: 10.1186/1471-2288-10-96)

set.seed(101)

### Load relevant packages
library(pROC)

### Define the C-statistic values for each model
### We will use the values from the internal validation, the reason for this is that then we can derive R2_CS_app, and adjust it for optimism
### in the calculation of N_min. Using the external validation would give us a direct estimate of R2_CS_adj, but it would only be valid in the population
### in which the model was validated, which was a chinese cohort. We assume we are developing the model in a British population.


### First define the number of events in each category (Benign = 1, Borderline malignent = 2, primary invaseive = 3, metastatic = 4)
EV1 <- 800
EV2 <- 55
EV3 <- 169
EV4 <- 42


### The C-statistic for Benign vs Borderline malignent is 0.82
### The C-statistic for Benign vs Prim invasive is 0.95
### The C-statistic for Benign vs Metastatic is 0.93
C.2 <- 0.82
C.3 <- 0.95
C.4 <- 0.93


### The outcome proportion for Benign vs Borderline malignent is 55/(800 + 55)
### The outcome proportion for Benign vs Prim invasive is 169/(800 + 169)
### The outcome proportion for Benign vs Metastatic is 42/(800 + 42)
Eprop.2 <- EV2/(EV1 + EV2)
Eprop.3 <- EV3/(EV1 + EV3)
Eprop.4 <- EV4/(EV1 + EV4)


### Create a function simulate a dataset as described in Riley et al, for a given C-statistic and event proportion
simulate.R2 <- function(N, prop.in, C.in){

  ## Create an empty dataset
  output.dat <- data.frame(matrix(ncol = 2, nrow = N))
  colnames(output.dat) <- c("Y", "LP")
  
  ## Create the outcome variable
  Y.vec <- rbinom(N, 1, prop.in)
  
  ## Create the vector of mean values for the linear predictor data generation
  Y.vec.mu <- Y.vec*sqrt(2)*qnorm(C.in, 0, 1)
  
  ## Generate the linear predictor
  LP.vec <- rnorm(N,Y.vec.mu,1)

  ## Assign these into an output dataset
  output.dat$Y <- as.integer(Y.vec)
  output.dat$LP <- LP.vec
  
  ## Fit a logistic regression to this dataset
  model.out <- glm(Y ~ LP.vec, family = binomial(link = "logit"), data = output.dat)
  model.out
  
  ## Check the AUC is correct, matches input
  #C.stat.sim <- as.numeric(roc(Y ~ LP.vec, data = output.dat)$auc)
  
  ## Fit a null model also, to calculate likelihood ratio
  model.null <- glm(Y ~ 1, family = binomial(link = "logit"), data = output.dat)
  
  ## Calculate likelihood ratio statistics
  LR <- as.numeric(-2*(logLik(model.null) - logLik(model.out)))

  ## Calculate R2_CS_APP
  R2_CS_APP <- 1 - exp(-LR/N)
  R2_CS_APP

  ## Output object
  return(R2_CS_APP)
}

### Calculate R2_CS according to the simulation approach of Riley.
R2_CS_sim.2 <- simulate.R2(1000000, Eprop.2, C.2)
R2_CS_sim.3 <- simulate.R2(1000000, Eprop.3, C.3)
R2_CS_sim.4 <- simulate.R2(1000000, Eprop.4, C.4)

R2_CS_sim.2
R2_CS_sim.3
R2_CS_sim.4

### Calculate N_DL_sim based on this assumed R2_CS value, to target a threshold of 0.9, by plugging into formula of Riley et al.
## Let p be number of predictors considered in initial model, supposing we will consider the same set before variable selection
p <- 16

## Let S be the level of shrinkage we are targeting
S <- 0.9

## Calculate N_DL_sim_temporary
N_DL_sim.2.temp <- p/((S - 1)*log(1 - R2_CS_sim.2/S))
N_DL_sim.3.temp <- p/((S - 1)*log(1 - R2_CS_sim.3/S))
N_DL_sim.4.temp <- p/((S - 1)*log(1 - R2_CS_sim.4/S))

N_DL_sim.2.temp
N_DL_sim.3.temp
N_DL_sim.4.temp

### Note that this is the number of individuals required for each distinct logistic regression model, but only a subset of the individuals
### recruited will be used in these models, so need to divide by the proportion of individuals that are used in each of the models
n.total <- EV1 + EV2 + EV3 + EV4
Eprop.total.2 <- (EV1 + EV2)/n.total
Eprop.total.3 <- (EV1 + EV3)/n.total
Eprop.total.4 <- (EV1 + EV4)/n.total

### Calculate N_DL_sim final
N_DL_sim.2 <- N_DL_sim.2.temp/Eprop.total.2
N_DL_sim.3 <- N_DL_sim.3.temp/Eprop.total.3
N_DL_sim.4 <- N_DL_sim.4.temp/Eprop.total.4

N_DL_sim.2.temp
N_DL_sim.3.temp
N_DL_sim.4.temp

N_DL_sim.2
N_DL_sim.3
N_DL_sim.4

