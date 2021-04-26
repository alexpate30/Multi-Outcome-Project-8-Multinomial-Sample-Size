### This program will summarise the results from the simulation
rm(list=ls())

## Load image
load("/mnt/bmh01-rds/vanstaa/practice2019/alex/MRC multistate/multinom simulation scen1to5 500000val.RData")


### Create lists of dataframes of the simulations output
scen1.dataframes <- vector("list", 9)
for (i in 1:9){scen1.dataframes[[i]] <- data.frame(scenario1.res[[(i)]])}

scen2.dataframes <- vector("list", 9)
for (i in 1:9){scen2.dataframes[[i]] <- data.frame(scenario2.res[[(i)]])}

scen3.dataframes <- vector("list", 9)
for (i in 1:9){scen3.dataframes[[i]] <- data.frame(scenario3.res[[(i)]])}

scen4.dataframes <- vector("list", 8)
for (i in 1:8){scen4.dataframes[[i]] <- data.frame(scenario4.res[[(i)]])}

scen5.dataframes <- vector("list", 8)
for (i in 1:8){scen5.dataframes[[i]] <- data.frame(scenario5.res[[(i)]])}


### Give names to list elements### Scenario 1
names(scen1.dataframes) <- c(100,250,500,1000,540,564,576,572,606)
names(scen2.dataframes) <- c(100,250,500,1000,569,699,628,570,571)
names(scen3.dataframes) <- c(100,250,500,1000,565,714,582,565,566)
names(scen4.dataframes) <- c(250,500,1000,647,1105,901,813,1025)
names(scen5.dataframes) <- c(250,500,1000,706,1742,1457,1156,2003)


sample.sizes1 <- c(100,250,500,1000,540,564,576,572,606)
sample.sizes2 <- c(100,250,500,1000,569,699,628,570,571)
sample.sizes3 <- c(100,250,500,1000,565,714,582,565,566)
sample.sizes4 <- c(250,500,1000,647,1105,901,813,1025)
sample.sizes5 <- c(250,500,1000,706,1742,1457,1156,2003)


## Want to calculate the heuristic shrinkage factor for each model
calc.S_VH <- function(data.in){
  data.out <- data.in
  for (i in 1:9){
      data.out[[i]]$S_VH_multinom <- 1 - 10/(data.out[[i]]$LR.multinom)
      data.out[[i]]$S_VH_seqlog1 <- 1 - 5/(data.out[[i]]$LR.seqlog.1)
      data.out[[i]]$S_VH_seqlog2 <- 1 - 5/(data.out[[i]]$LR.seqlog.2)
      data.out[[i]]$S_VH_dislog1 <- 1 - 5/(data.out[[i]]$LR.dislog.1)
      data.out[[i]]$S_VH_dislog2 <- 1 - 5/(data.out[[i]]$LR.dislog.2)
    }
  return(data.out)
}

calc.S_VH.i8 <- function(data.in){
  data.out <- data.in
  for (i in 1:8){
    data.out[[i]]$S_VH_multinom <- 1 - 10/(data.out[[i]]$LR.multinom)
    data.out[[i]]$S_VH_seqlog1 <- 1 - 5/(data.out[[i]]$LR.seqlog.1)
    data.out[[i]]$S_VH_seqlog2 <- 1 - 5/(data.out[[i]]$LR.seqlog.2)
    data.out[[i]]$S_VH_dislog1 <- 1 - 5/(data.out[[i]]$LR.dislog.1)
    data.out[[i]]$S_VH_dislog2 <- 1 - 5/(data.out[[i]]$LR.dislog.2)
  }
  return(data.out)
}


## Apply functions
scen1.dataframes.S <- calc.S_VH(scen1.dataframes)
scen2.dataframes.S <- calc.S_VH(scen2.dataframes)
scen3.dataframes.S <- calc.S_VH(scen3.dataframes)
scen4.dataframes.S <- calc.S_VH.i8(scen4.dataframes)
scen5.dataframes.S <- calc.S_VH.i8(scen5.dataframes)




### Calculate R2_CS_APP and R2_CS_ADJ for each model
### Going to just write a separate function for each as they all have a different set of sample sizes
for (i in 1:9){scen1.dataframes.S[[i]]$R2_CS_APP <- 1 - exp(-scen1.dataframes.S[[i]]$LR.multinom/sample.sizes1[i])
               scen1.dataframes.S[[i]]$R2_CS_ADJ <- scen1.dataframes.S[[i]]$R2_CS_APP*scen1.dataframes.S[[i]]$S_VH_multinom}
for (i in 1:9){scen2.dataframes.S[[i]]$R2_CS_APP <- 1 - exp(-scen2.dataframes.S[[i]]$LR.multinom/sample.sizes2[i])
               scen2.dataframes.S[[i]]$R2_CS_ADJ <- scen2.dataframes.S[[i]]$R2_CS_APP*scen2.dataframes.S[[i]]$S_VH_multinom}
for (i in 1:9){scen3.dataframes.S[[i]]$R2_CS_APP <- 1 - exp(-scen3.dataframes.S[[i]]$LR.multinom/sample.sizes3[i])
               scen3.dataframes.S[[i]]$R2_CS_ADJ <- scen3.dataframes.S[[i]]$R2_CS_APP*scen3.dataframes.S[[i]]$S_VH_multinom}
for (i in 1:8){scen4.dataframes.S[[i]]$R2_CS_APP <- 1 - exp(-scen4.dataframes.S[[i]]$LR.multinom/sample.sizes4[i])
               scen4.dataframes.S[[i]]$R2_CS_ADJ <- scen4.dataframes.S[[i]]$R2_CS_APP*scen4.dataframes.S[[i]]$S_VH_multinom}
for (i in 1:8){scen5.dataframes.S[[i]]$R2_CS_APP <- 1 - exp(-scen5.dataframes.S[[i]]$LR.multinom/sample.sizes5[i])
               scen5.dataframes.S[[i]]$R2_CS_ADJ <- scen5.dataframes.S[[i]]$R2_CS_APP*scen5.dataframes.S[[i]]$S_VH_multinom}


### Calculate median of each entity of interest and put into a table
get.output.tables.med <- function(data.in){
  tables.out <- rbind(
    # First row
    c(capture.output(cat(round(quantile(data.in[[1]]$Slope.devel.m.cali.m.1, probs = c(0.5)),3), " (", round(mean(data.in[[1]]$Slope.devel.m.cali.m.1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[1]]$Slope.devel.m.cali.m.2, probs = c(0.5)),3), " (", round(mean(data.in[[1]]$Slope.devel.m.cali.m.2),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[1]]$S_VH_multinom, probs = c(0.5)),3), " (", round(mean(data.in[[1]]$S_VH_multinom),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[1]]$Slope.devel.d.cali.d.1, probs = c(0.5)),3), " (", round(mean(data.in[[1]]$Slope.devel.d.cali.d.1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[1]]$Slope.devel.d.cali.d.2, probs = c(0.5)),3), " (", round(mean(data.in[[1]]$Slope.devel.d.cali.d.2),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[1]]$S_VH_dislog1, probs = c(0.5)),3), " (", round(mean(data.in[[1]]$S_VH_dislog1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[1]]$S_VH_dislog2, probs = c(0.5)),3), " (", round(mean(data.in[[1]]$S_VH_dislog2),3), ")", sep = ""))),
    # Second row
    c(capture.output(cat(round(quantile(data.in[[2]]$Slope.devel.m.cali.m.1, probs = c(0.5)),3), " (", round(mean(data.in[[2]]$Slope.devel.m.cali.m.1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[2]]$Slope.devel.m.cali.m.2, probs = c(0.5)),3), " (", round(mean(data.in[[2]]$Slope.devel.m.cali.m.2),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[2]]$S_VH_multinom, probs = c(0.5)),3), " (", round(mean(data.in[[2]]$S_VH_multinom),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[2]]$Slope.devel.d.cali.d.1, probs = c(0.5)),3), " (", round(mean(data.in[[2]]$Slope.devel.d.cali.d.1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[2]]$Slope.devel.d.cali.d.2, probs = c(0.5)),3), " (", round(mean(data.in[[2]]$Slope.devel.d.cali.d.2),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[2]]$S_VH_dislog1, probs = c(0.5)),3), " (", round(mean(data.in[[2]]$S_VH_dislog1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[2]]$S_VH_dislog2, probs = c(0.5)),3), " (", round(mean(data.in[[2]]$S_VH_dislog2),3), ")", sep = ""))),
    # Third row
    c(capture.output(cat(round(quantile(data.in[[3]]$Slope.devel.m.cali.m.1, probs = c(0.5)),3), " (", round(mean(data.in[[3]]$Slope.devel.m.cali.m.1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[3]]$Slope.devel.m.cali.m.2, probs = c(0.5)),3), " (", round(mean(data.in[[3]]$Slope.devel.m.cali.m.2),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[3]]$S_VH_multinom, probs = c(0.5)),3), " (", round(mean(data.in[[3]]$S_VH_multinom),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[3]]$Slope.devel.d.cali.d.1, probs = c(0.5)),3), " (", round(mean(data.in[[3]]$Slope.devel.d.cali.d.1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[3]]$Slope.devel.d.cali.d.2, probs = c(0.5)),3), " (", round(mean(data.in[[3]]$Slope.devel.d.cali.d.2),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[3]]$S_VH_dislog1, probs = c(0.5)),3), " (", round(mean(data.in[[3]]$S_VH_dislog1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[3]]$S_VH_dislog2, probs = c(0.5)),3), " (", round(mean(data.in[[3]]$S_VH_dislog2),3), ")", sep = ""))),
    # Fourth row
    c(capture.output(cat(round(quantile(data.in[[4]]$Slope.devel.m.cali.m.1, probs = c(0.5)),3), " (", round(mean(data.in[[4]]$Slope.devel.m.cali.m.1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[4]]$Slope.devel.m.cali.m.2, probs = c(0.5)),3), " (", round(mean(data.in[[4]]$Slope.devel.m.cali.m.2),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[4]]$S_VH_multinom, probs = c(0.5)),3), " (", round(mean(data.in[[4]]$S_VH_multinom),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[4]]$Slope.devel.d.cali.d.1, probs = c(0.5)),3), " (", round(mean(data.in[[4]]$Slope.devel.d.cali.d.1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[4]]$Slope.devel.d.cali.d.2, probs = c(0.5)),3), " (", round(mean(data.in[[4]]$Slope.devel.d.cali.d.2),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[4]]$S_VH_dislog1, probs = c(0.5)),3), " (", round(mean(data.in[[4]]$S_VH_dislog1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[4]]$S_VH_dislog2, probs = c(0.5)),3), " (", round(mean(data.in[[4]]$S_VH_dislog2),3), ")", sep = ""))),
    # Fifth row
    c(capture.output(cat(round(quantile(data.in[[5]]$Slope.devel.m.cali.m.1, probs = c(0.5)),3), " (", round(mean(data.in[[5]]$Slope.devel.m.cali.m.1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[5]]$Slope.devel.m.cali.m.2, probs = c(0.5)),3), " (", round(mean(data.in[[5]]$Slope.devel.m.cali.m.2),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[5]]$S_VH_multinom, probs = c(0.5)),3), " (", round(mean(data.in[[5]]$S_VH_multinom),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[5]]$Slope.devel.d.cali.d.1, probs = c(0.5)),3), " (", round(mean(data.in[[5]]$Slope.devel.d.cali.d.1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[5]]$Slope.devel.d.cali.d.2, probs = c(0.5)),3), " (", round(mean(data.in[[5]]$Slope.devel.d.cali.d.2),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[5]]$S_VH_dislog1, probs = c(0.5)),3), " (", round(mean(data.in[[5]]$S_VH_dislog1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[5]]$S_VH_dislog2, probs = c(0.5)),3), " (", round(mean(data.in[[5]]$S_VH_dislog2),3), ")", sep = ""))),
    # Sixth row
    c(capture.output(cat(round(quantile(data.in[[6]]$Slope.devel.m.cali.m.1, probs = c(0.5)),3), " (", round(mean(data.in[[6]]$Slope.devel.m.cali.m.1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[6]]$Slope.devel.m.cali.m.2, probs = c(0.5)),3), " (", round(mean(data.in[[6]]$Slope.devel.m.cali.m.2),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[6]]$S_VH_multinom, probs = c(0.5)),3), " (", round(mean(data.in[[6]]$S_VH_multinom),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[6]]$Slope.devel.d.cali.d.1, probs = c(0.5)),3), " (", round(mean(data.in[[6]]$Slope.devel.d.cali.d.1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[6]]$Slope.devel.d.cali.d.2, probs = c(0.5)),3), " (", round(mean(data.in[[6]]$Slope.devel.d.cali.d.2),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[6]]$S_VH_dislog1, probs = c(0.5)),3), " (", round(mean(data.in[[6]]$S_VH_dislog1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[6]]$S_VH_dislog2, probs = c(0.5)),3), " (", round(mean(data.in[[6]]$S_VH_dislog2),3), ")", sep = ""))),
    # Seventh row
    c(capture.output(cat(round(quantile(data.in[[7]]$Slope.devel.m.cali.m.1, probs = c(0.5)),3), " (", round(mean(data.in[[7]]$Slope.devel.m.cali.m.1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[7]]$Slope.devel.m.cali.m.2, probs = c(0.5)),3), " (", round(mean(data.in[[7]]$Slope.devel.m.cali.m.2),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[7]]$S_VH_multinom, probs = c(0.5)),3), " (", round(mean(data.in[[7]]$S_VH_multinom),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[7]]$Slope.devel.d.cali.d.1, probs = c(0.5)),3), " (", round(mean(data.in[[7]]$Slope.devel.d.cali.d.1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[7]]$Slope.devel.d.cali.d.2, probs = c(0.5)),3), " (", round(mean(data.in[[7]]$Slope.devel.d.cali.d.2),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[7]]$S_VH_dislog1, probs = c(0.5)),3), " (", round(mean(data.in[[7]]$S_VH_dislog1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[7]]$S_VH_dislog2, probs = c(0.5)),3), " (", round(mean(data.in[[7]]$S_VH_dislog2),3), ")", sep = "")))
    )
    return(tables.out)}
      

### Repeat but for the tables that just have 6 rows
get.output.tables.med.i6 <- function(data.in){
  tables.out <- rbind(
    # First row
    c(capture.output(cat(round(quantile(data.in[[1]]$Slope.devel.m.cali.m.1, probs = c(0.5)),3), " (", round(mean(data.in[[1]]$Slope.devel.m.cali.m.1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[1]]$Slope.devel.m.cali.m.2, probs = c(0.5)),3), " (", round(mean(data.in[[1]]$Slope.devel.m.cali.m.2),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[1]]$S_VH_multinom, probs = c(0.5)),3), " (", round(mean(data.in[[1]]$S_VH_multinom),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[1]]$Slope.devel.d.cali.d.1, probs = c(0.5)),3), " (", round(mean(data.in[[1]]$Slope.devel.d.cali.d.1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[1]]$Slope.devel.d.cali.d.2, probs = c(0.5)),3), " (", round(mean(data.in[[1]]$Slope.devel.d.cali.d.2),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[1]]$S_VH_dislog1, probs = c(0.5)),3), " (", round(mean(data.in[[1]]$S_VH_dislog1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[1]]$S_VH_dislog2, probs = c(0.5)),3), " (", round(mean(data.in[[1]]$S_VH_dislog2),3), ")", sep = ""))),
    # Second row
    c(capture.output(cat(round(quantile(data.in[[2]]$Slope.devel.m.cali.m.1, probs = c(0.5)),3), " (", round(mean(data.in[[2]]$Slope.devel.m.cali.m.1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[2]]$Slope.devel.m.cali.m.2, probs = c(0.5)),3), " (", round(mean(data.in[[2]]$Slope.devel.m.cali.m.2),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[2]]$S_VH_multinom, probs = c(0.5)),3), " (", round(mean(data.in[[2]]$S_VH_multinom),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[2]]$Slope.devel.d.cali.d.1, probs = c(0.5)),3), " (", round(mean(data.in[[2]]$Slope.devel.d.cali.d.1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[2]]$Slope.devel.d.cali.d.2, probs = c(0.5)),3), " (", round(mean(data.in[[2]]$Slope.devel.d.cali.d.2),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[2]]$S_VH_dislog1, probs = c(0.5)),3), " (", round(mean(data.in[[2]]$S_VH_dislog1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[2]]$S_VH_dislog2, probs = c(0.5)),3), " (", round(mean(data.in[[2]]$S_VH_dislog2),3), ")", sep = ""))),
    # Third row
    c(capture.output(cat(round(quantile(data.in[[3]]$Slope.devel.m.cali.m.1, probs = c(0.5)),3), " (", round(mean(data.in[[3]]$Slope.devel.m.cali.m.1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[3]]$Slope.devel.m.cali.m.2, probs = c(0.5)),3), " (", round(mean(data.in[[3]]$Slope.devel.m.cali.m.2),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[3]]$S_VH_multinom, probs = c(0.5)),3), " (", round(mean(data.in[[3]]$S_VH_multinom),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[3]]$Slope.devel.d.cali.d.1, probs = c(0.5)),3), " (", round(mean(data.in[[3]]$Slope.devel.d.cali.d.1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[3]]$Slope.devel.d.cali.d.2, probs = c(0.5)),3), " (", round(mean(data.in[[3]]$Slope.devel.d.cali.d.2),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[3]]$S_VH_dislog1, probs = c(0.5)),3), " (", round(mean(data.in[[3]]$S_VH_dislog1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[3]]$S_VH_dislog2, probs = c(0.5)),3), " (", round(mean(data.in[[3]]$S_VH_dislog2),3), ")", sep = ""))),
    # Fourth row
    c(capture.output(cat(round(quantile(data.in[[4]]$Slope.devel.m.cali.m.1, probs = c(0.5)),3), " (", round(mean(data.in[[4]]$Slope.devel.m.cali.m.1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[4]]$Slope.devel.m.cali.m.2, probs = c(0.5)),3), " (", round(mean(data.in[[4]]$Slope.devel.m.cali.m.2),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[4]]$S_VH_multinom, probs = c(0.5)),3), " (", round(mean(data.in[[4]]$S_VH_multinom),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[4]]$Slope.devel.d.cali.d.1, probs = c(0.5)),3), " (", round(mean(data.in[[4]]$Slope.devel.d.cali.d.1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[4]]$Slope.devel.d.cali.d.2, probs = c(0.5)),3), " (", round(mean(data.in[[4]]$Slope.devel.d.cali.d.2),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[4]]$S_VH_dislog1, probs = c(0.5)),3), " (", round(mean(data.in[[4]]$S_VH_dislog1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[4]]$S_VH_dislog2, probs = c(0.5)),3), " (", round(mean(data.in[[4]]$S_VH_dislog2),3), ")", sep = ""))),
    # Fifth row
    c(capture.output(cat(round(quantile(data.in[[5]]$Slope.devel.m.cali.m.1, probs = c(0.5)),3), " (", round(mean(data.in[[5]]$Slope.devel.m.cali.m.1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[5]]$Slope.devel.m.cali.m.2, probs = c(0.5)),3), " (", round(mean(data.in[[5]]$Slope.devel.m.cali.m.2),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[5]]$S_VH_multinom, probs = c(0.5)),3), " (", round(mean(data.in[[5]]$S_VH_multinom),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[5]]$Slope.devel.d.cali.d.1, probs = c(0.5)),3), " (", round(mean(data.in[[5]]$Slope.devel.d.cali.d.1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[5]]$Slope.devel.d.cali.d.2, probs = c(0.5)),3), " (", round(mean(data.in[[5]]$Slope.devel.d.cali.d.2),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[5]]$S_VH_dislog1, probs = c(0.5)),3), " (", round(mean(data.in[[5]]$S_VH_dislog1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[5]]$S_VH_dislog2, probs = c(0.5)),3), " (", round(mean(data.in[[5]]$S_VH_dislog2),3), ")", sep = ""))),
    # Sixth row
    c(capture.output(cat(round(quantile(data.in[[6]]$Slope.devel.m.cali.m.1, probs = c(0.5)),3), " (", round(mean(data.in[[6]]$Slope.devel.m.cali.m.1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[6]]$Slope.devel.m.cali.m.2, probs = c(0.5)),3), " (", round(mean(data.in[[6]]$Slope.devel.m.cali.m.2),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[6]]$S_VH_multinom, probs = c(0.5)),3), " (", round(mean(data.in[[6]]$S_VH_multinom),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[6]]$Slope.devel.d.cali.d.1, probs = c(0.5)),3), " (", round(mean(data.in[[6]]$Slope.devel.d.cali.d.1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[6]]$Slope.devel.d.cali.d.2, probs = c(0.5)),3), " (", round(mean(data.in[[6]]$Slope.devel.d.cali.d.2),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[6]]$S_VH_dislog1, probs = c(0.5)),3), " (", round(mean(data.in[[6]]$S_VH_dislog1),3), ")", sep = "")),
      capture.output(cat(round(quantile(data.in[[6]]$S_VH_dislog2, probs = c(0.5)),3), " (", round(mean(data.in[[6]]$S_VH_dislog2),3), ")", sep = "")))
  )
  return(tables.out)}

### And apply it to each set of results
output.scen1.med <- get.output.tables.med(scen1.dataframes.S)
output.scen2.med <- get.output.tables.med(scen2.dataframes.S)
output.scen3.med <- get.output.tables.med(scen3.dataframes.S)
output.scen4.med <- get.output.tables.med.i6(scen4.dataframes.S)
output.scen5.med <- get.output.tables.med.i6(scen5.dataframes.S)


save.image("R_out/multinom results Table2.RData")
