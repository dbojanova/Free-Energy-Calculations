library(CHNOSZ)
library(zoo)
library(tidyverse)

free_energy_calcs <- function(redox, reactions, reaction_names, minerals, I) {
  
  # ~~~ Extract temperature (C) and pressure (bar) values ~~~
  t <- redox$Temperature
  p <- redox$Pressure
  pH <- redox$pH
  redox$`H+` <- 10^(-redox$pH)
  
  redox_spp <- names(redox)[5:length(redox)]
  
  E.units("J")
  
  Gknot <- list()
  
  for(i in 1:length(reactions)){
    state <- rep("aq",length(reactions[[i]][[1]]))
    for(x in 1:length(reactions[[i]][[1]])){
      if(any(minerals %in% reactions[[i]][[1]][x])){
        state[x] <- 'cr'
      }
    }
    Gknot[[i]] <- subcrt(reactions[[i]][[1]],state,reactions[[i]][[2]],T=t, P=p)$out$G/1000
  }
  
  ionic = data.frame(
    ionicStrength = c(0.001,0.001,0.001,0.001,0.001,0.01,0.01,0.01,0.01,0.01,0.1,0.1,0.1,0.1,0.1,0.7,0.7,0.7,0.7,0.7),
    temperature = c(0,25,50,75,100,0,25,50,75,100,0,25,50,75,100,0,25,50,75,100),
    minus3 = c(0.73,0.72,0.71,0.70,0.69,0.4,0.39,0.38,0.36,0.34,0.1,0.09,0.09,0.08,0.07,0.02,0.02,0.01,0.01,0.01),
    minus2 = c(0.87,0.87,0.86,0.85,0.85,0.67,0.66,0.65,0.63,0.62,0.36,0.35,0.34,0.32,0.30,0.17,0.16,0.15,0.14,0.12),
    minus1 = c(0.97,0.96,0.96,0.96,0.96,0.90,0.90,0.90,0.89,0.89,0.78,0.77,0.77,0.76,0.74,0.67,0.66,0.65,0.64,0.62),
    zero = c(rep(1,15),rep(0.99,5)),
    plus1 = c(0.97,0.96,0.96,0.96,0.96,0.90,0.90,0.90,0.89,0.89,0.78,0.77,0.77,0.76,0.74,0.67,0.66,0.65,0.64,0.62),
    plus2 = c(0.87,0.87,0.86,0.86,0.85,0.68,0.68,0.66,0.65,0.63,0.41,0.40,0.39,0.37,0.35,0.25,0.24,0.23,0.21,0.19),
    plus3 = c(0.74,0.74,0.73,0.71,0.70,0.45,0.44,0.43,0.41,0.39,0.19,0.18,0.17,0.15,0.14,0.09,0.08,0.08,0.07,0.06)
  )
  
  ionicUsed = subset(ionic, ionic$ionicStrength == I)
  
  intercept = c()
  slope = c()
  
  for(temp in 3:9){
    fit = lm(formula = as.matrix(ionicUsed[temp])~ionicUsed$temperature)
    intercept[temp - 2] = fit$coefficients[1]
    slope[temp - 2] = fit$coefficients[2]
  }
  
  for(conc in 5:length(redox)){
    lastChar = substring(redox_spp[conc-4], nchar(redox_spp[conc-4]))
    lastTwoChar = substring(redox_spp[conc-4], nchar(redox_spp[conc-4])-1,nchar(redox_spp[conc-4]))
    
    if(lastTwoChar == "-3"){
      redox[conc] = redox[conc]*(slope[1]*t + intercept[1])
    }else if(lastTwoChar == "-2"){
      redox[conc] = redox[conc]*(slope[2]*t + intercept[2])
    }else if(lastChar == "-"){
      redox[conc] = redox[conc]*(slope[3]*t + intercept[3])
    }else if(lastChar == "+"){
      redox[conc] = redox[conc]*(slope[5]*t + intercept[5])
    }else if(lastTwoChar == "+2"){
      redox[conc] = redox[conc]*(slope[6]*t + intercept[6])
    }else if(lastTwoChar == "+3"){
      redox[conc] = redox[conc]*(slope[7]*t + intercept[7])
    }else{
      redox[conc] = redox[conc]*(slope[4]*t + intercept[4])
    }
  }
  
  activity <- list()
  
  for(c in 1:length(reactions)){
    spp_placement <- list()
    for(d in 1:length(redox_spp)){
      v <- reactions[[c]][[1]]
      if(redox_spp[d] %in% v){
        spp_placement[[match(redox_spp[d], v)]] <- redox[which(names(redox)==redox_spp[d])]
      }
    }
    activity[[c]] <- spp_placement
  }
  
  multiplyList <- function(myList){
    m <- 1
    for(list in 1:length(myList)){
      m <- m * myList[[list]]
    }
    return(m)
  }
  
  Q <- list()
  
  for(r in 1:length(reactions)){
    A_react <- list()
    for(num in 1:length(reactions[[r]][[2]])){
      A_react[[num]] <- activity[[r]][[num]] ** reactions[[r]][[2]][[num]]
    }
    Q[[r]] <- multiplyList(A_react)
  }
  
  RealGibbs <- data.frame("Sample" = redox$Sample, "Temp" = t, "pH" = pH)
  
  for(gibbs in 1:length(reactions)){
    RealGibbs[gibbs+3] <- ((Gknot[[gibbs]] * Q[[gibbs]]) / Q[[gibbs]]) + 0.0083145 * (t + 273.15) * log(Q[[gibbs]])
    names(RealGibbs)[gibbs+3] <- reaction_names[gibbs]
  }
  
  return(RealGibbs)
}
