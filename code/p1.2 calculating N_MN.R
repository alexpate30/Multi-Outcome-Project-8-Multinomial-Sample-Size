### This program will run a small simulation where I simulate multinomial logistic data and fit models, and assess the calibration slopes for various sample sizes
#install.packages("VGAM")

library(VGAM)

### Set seed
set.seed(101)

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


### Before running the simulation, I will create a validation dataset of size 100,000 
### for each scenario (set of coeffs)
data.valid.s1 <- create.dataset(500000, 0, 1, -0.5, -0.25, 0.5, 0.75,
                             0, 0.75, -1, -0.5, -0.75, 0.25)
data.valid.s2 <- create.dataset(500000, 0, 1, -0.5, -0.25, 0.5, 0.75,
                                -1, 0.75, -1, -0.5, -0.75, 0.25)
data.valid.s3 <- create.dataset(500000, -0.4, 1, -0.5, -0.25, 0.5, 0.75,
                                -1, 0.75, -1, -0.5, -0.75, 0.25)
data.valid.s4 <- create.dataset(500000, -0.5, 1, -0.5, -0.25, 0.5, 0.75,
                                -2, 0.75, -1, -0.5, -0.75, 0.25)


### Now fit a model in each of these and get the loglikelihood
model.s1 <- vgam(y.cat ~ x1 + x2 + x3 + x4 + x5, family = multinomial(ref = "cat1"), data = data.valid.s1)
model.s2 <- vgam(y.cat ~ x1 + x2 + x3 + x4 + x5, family = multinomial(ref = "cat1"), data = data.valid.s2)
model.s3 <- vgam(y.cat ~ x1 + x2 + x3 + x4 + x5, family = multinomial(ref = "cat1"), data = data.valid.s3)
model.s4 <- vgam(y.cat ~ x1 + x2 + x3 + x4 + x5, family = multinomial(ref = "cat1"), data = data.valid.s4)

model.null.s1 <- vgam(y.cat ~ 1, family = multinomial(ref = "cat1"), data = data.valid.s1)
model.null.s2 <- vgam(y.cat ~ 1, family = multinomial(ref = "cat1"), data = data.valid.s2)
model.null.s3 <- vgam(y.cat ~ 1, family = multinomial(ref = "cat1"), data = data.valid.s3)
model.null.s4 <- vgam(y.cat ~ 1, family = multinomial(ref = "cat1"), data = data.valid.s4)

### Extract  loglikelihoods
lr.1 <- -2*(model.null.s1@criterion$loglikelihood - model.s1@criterion$loglikelihood)
lr.2 <- -2*(model.null.s2@criterion$loglikelihood - model.s2@criterion$loglikelihood)
lr.3 <- -2*(model.null.s3@criterion$loglikelihood - model.s3@criterion$loglikelihood)
lr.4 <- -2*(model.null.s4@criterion$loglikelihood - model.s4@criterion$loglikelihood)

### Calculaet R2_APP
R2_APP.1 <- 1 - exp(-lr.1/500000)
R2_APP.2 <- 1 - exp(-lr.2/500000)
R2_APP.3 <- 1 - exp(-lr.3/500000)
R2_APP.4 <- 1 - exp(-lr.4/500000)

### Calculate S_VH
S_VH.1 <- 1 + (10/(500000*log(1-R2_APP.1)))
S_VH.2 <- 1 + (10/(500000*log(1-R2_APP.2)))
S_VH.3 <- 1 + (10/(500000*log(1-R2_APP.3)))
S_VH.4 <- 1 + (10/(500000*log(1-R2_APP.4)))

### Calculate R2_ADJ
R2_ADJ.1 <- S_VH.1*R2_APP.1
R2_ADJ.2 <- S_VH.1*R2_APP.2
R2_ADJ.3 <- S_VH.1*R2_APP.3
R2_ADJ.4 <- S_VH.1*R2_APP.4


### Calculate required n
n.req.1 <- 10/((0.9-1)*log(1-(R2_ADJ.1/0.9)))
n.req.2 <- 10/((0.9-1)*log(1-(R2_ADJ.2/0.9)))
n.req.3 <- 10/((0.9-1)*log(1-(R2_ADJ.3/0.9)))
n.req.4 <- 10/((0.9-1)*log(1-(R2_ADJ.4/0.9)))
      
### Print n.req
n.req.1
n.req.2
n.req.3
n.req.4
          
rm(data.valid.s1, data.valid.s2, data.valid.s3, data.valid.s4, model.s1, model.s2, model.s3, model.s4,
   model.null.s1, model.null.s2, model.null.s3, model.null.s4)

save.image("R_out/calculating N_MN.RData")
