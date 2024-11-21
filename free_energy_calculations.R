#Free energy calculations (real delta G) 

library(CHNOSZ)
library(plyr)
library(zoo)

#set the working directory (where you have saved the physicochemical_summary.csv)

  setwd("")

#import site redox - includes temperature, pressure, and concentrations (M)

  redox <- read.csv("physicochemical_summary.csv",check.names=FALSE)
  
#remove rows with NA for depth and methane
  redox <- redox[rowSums(is.na(redox[,4:5])) == 0, ]  

#replace NA values with previous value
  redox$Pressure <- na.locf(redox$Pressure)

#pull out the redox species names and indicate minerals of interest
#MUST write minerals, if using any
  names(redox) <- c("Site","Hole","Core","Depth","Temperature","Pressure","HCO3-",
                    "Br-","Li+2","H+","P","salinity","SO4-2","S",
                    "NH4+","PO4-3","HS-","Ba+2","B","Ca+2","Cl-",
                    "Fe+2","Mg+2","Mn+2","K+","Si+4","Na+","Sr+2",
                    "CO2","H2O","MnO2","FeOOH","Fe(OH)3","Fe2O3","Fe3O4",
                    "CH4","ethane","propane","butane","pentane","hexane","H2","CO")
  
#isolate redox species names and minerals of interest
  
  redox_spp <- names(redox)[7:length(redox)]

  
  minerals <- c("MnO2","FeOOH","Fe(OH)3","Fe2O3","Fe3O4")
  
#Select the ionic strength of your environment (either 0.001, 0,01, 0.1, or 0.7)
  
  I = 0.7
  
#Set up reactions
  # Reaction names can be whatever you wish, just keep them consistent
  # Format is: 
    #reactionName = list(c(reactant1, reactant 2, product1, product2),c(-coef1, -coef2, coef3, coef4)
    #the negative coefficient indicates that the respective species is a reactant
  
  #methanogenesis
    Hmethano <- list(c("H2","CO2","CH4","H2O"),c(-4,-1,1,2))
    HMmethano <- list(c("H2","methanol","CH4","H2O"),c(-1,-1,1,1))
    Amethano <- list(c("H+","CH3COO-","CH4","CO2"),c(-1,-1,1,1))
    Fmethano <- list(c("H+","HCO2-","CH4","CO2","H2O"),c(-4,-4,1,3,2))
    Mmethano <- list(c("methanol","CH4","CO2","H2O"),c(-4,3,1,2))
    methyl_methano <- list(c("methanamine","H2O","H+","CH4","CO2","NH4+"),
                           c(-4,-2,-4,3,1,4))
    Dmethyl_methano <- list(c("dimethylamine","H2O","H+","CH4","CO2","NH4+"),
                            c(-2,-2,-2,3,1,2))
    Tmethyl_methano <- list(c("trimethylamine","H2O","H+","CH4","CO2","NH4+"),
                            c(-4,-6,-4,9,3,4))
    DMSmethano <- list(c("dimethyl sulfide","H2O","CH4","CO2","HS-","H+"),
                       c(-2,-2,3,1,2,2))
    
  reactions <- list(Hmethano,Amethano,Fmethano,HMmethano,Mmethano,methyl_methano
                    ,Dmethyl_methano,Tmethyl_methano,DMSmethano)

  reaction_names <- c("Hmethano","Amethano","Fmethano","HMmethano","Mmethano",
                      "methyl_methano","Dmethyl_methano","Tmethyl_methano",
                      "DMSmethano")
  

########################################################################################################
                               ### WORKING CODE - DO NOT CHANGE ###
########################################################################################################

#Extract temperature (C) and pressure (bar) values
  
  t <- redox$Temperature
  p <- redox$Pressure

#Calculate delta G knot
  
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
 

#Calculate activites using activity coeffient estimates for species charge from Amend and LaRowe (2019)
  
  #Summary of activity coeffiencients for each ionic strength
  
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
  
  #Linear fit for ionic strength of choice and temperatures
  
  ionicUsed = subset(ionic, ionic$ionicStrength == I) 
  
  intercept = c()
  slope = c()
  
  for(temp in 3:9){
    fit = lm(formula = as.matrix(ionicUsed[temp])~ionicUsed$temperature)
    intercept[temp - 2] = fit$coefficients[1]
    slope[temp - 2] = fit$coefficients[2]
  }

  #Calculations based on the species in the original input redox

  for(conc in 7:length(redox)){
    lastChar = substring(redox_spp[conc-6], nchar(redox_spp[conc-6]))
    lastTwoChar = substring(redox_spp[conc-6], nchar(redox_spp[conc-6])-1,nchar(redox_spp[conc-6]))
    
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
  head(redox)
#Assign activities
  
  activity <- list()
  
  for(c in 1:length(reactions)){
    spp_placement <- list()
    for(d in 1:length(redox_spp)){
      v<-reactions[[c]][[1]]
      if(redox_spp[d] %in% v){
        spp_placement[[match(redox_spp[d], v)]] <- redox[which(names(redox)==redox_spp[d])]
      }
    }
    activity[[c]]<- spp_placement
  }

#Make multiplication function for lists
  
  multiplyList <- function(myList){
    m <- 1
    for(list in 1:length(myList)){
      m <- m*myList[[list]]
    }
    return(m)
  }
  
#Calculate Q  
  
Q <- list()

for(r in 1:length(reactions)){
  A_react <- list()
  for(num in 1:length(reactions[[r]][[2]])){
    A_react[[num]] <- activity[[r]][[num]]**reactions[[r]][[2]][[num]]
  }
  Q[[r]] <- multiplyList(A_react)
}

#Calculate the real Gibbs Energy for every reaction and return in redoxframe

RealGibbs <- data.frame("Depth" =redox$Depth,"Temp" = t,"Pressure" = p)


for(gibbs in 1:length(reactions)){
  RealGibbs[gibbs+3] <- ((Gknot[[gibbs]]*Q[[gibbs]])/Q[[gibbs]])+0.0083145*(t+273.15)*log(Q[[gibbs]])
  names(RealGibbs)[gibbs+3] <- reaction_names[gibbs]
}

#Export 
write.csv(RealGibbs, file = "Guaymas_RealGibbs_methano.csv",row.names = FALSE) 


