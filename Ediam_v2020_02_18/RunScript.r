#in pc
  root<-"~/Documents/macroeconomia/proyectoFinal"
  model.version<-"Ediam_v2020_02_18"
  dir.calib <- file.path(root, model.version, "CalibrationScripts")
  dir.model <- file.path(root, model.version)
  dir.data <- file.path(root, model.version, "Data")

####

#Load Data with scenarios

 #Select Population Scenario
#  pop.scenario.name<-"UN.Median.PI"
  pop.scenario.name<-"UN.Lower95.PI"
 #Load population scenario
  pop.scenario<<-read.csv(file.path(dir.data, "PopulationScenarios.csv"))
  pop.scenario<<-subset(pop.scenario,pop.scenario$Scenario==pop.scenario.name)
 #Select price of oil scenario
  oil.scenario.name<-"Baseline"
 #Load oil scenario
  oil.scenario<<-read.csv(file.path(dir.data, "OilPriceScenarios.csv"))
  oil.scenario<<-subset(oil.scenario,oil.scenario$Scenario==oil.scenario.name)
#Select climate change scenario
  climate.scenario<<-read.csv(file.path(dir.data, "ClimateScenarios.csv"))
  #climate.scenario.name<-"bcc-csm1-1"
  climate.scenario.name<-"GFDL-ESM2G"
  climate.scenario<<-subset(climate.scenario,climate.scenario$Climate.Model==climate.scenario.name)

#Load script for running Ediam
  source(file.path(dir.model, "ediam_Main.r"))
##Load script for plotting Ediam results
  source(file.path(dir.model, "ediam_Plot.r"))
#Load script for determining initial conditions
  source(file.path(dir.calib, "ediam_InitialConditions.r"))
#Load EDIAM Equations
  source(file.path(dir.model, "ediam_Equations_discrete.r"))

# Specify function for evaluating the performance of policy interventi
#definitely a maximization problem

ediamPolicyEval<-function(x,verbose=FALSE){

#Simulation Time Step
  TimeStep<<-5
#Base year for simulation: 2014
#historic data with hp
  Y_re.Nh <<- 490.0097 #Mtoe
  Y_ce.Nh <<- 4255.120 #Mtoe
  Y_re.Sh <<- 1386.1706 #Mtoe
  Y_ce.Sh <<- 7002.889 #Mtoe

#Oil Supply
  Re.Nh<<- 4255.120 #Mtoe
  Re.Sh<<- 7002.889 #Mtoe
  ReToGDP.Nh<<- 0.06227680 # [1]
  ReToGDP.Sh<<- 0.18765958 # [1]

#Oil prices
   oil.times<<- oil.scenario$Year - 2014
   Price.oil.y<<- oil.scenario$Price.Oil

#GDP with hp
    GDP.Nh<<- 47034.07 #(billion 2010 USD using exchange rates)
    GDP.Sh<<- 25688.192 #(billion 2010 USD using exchange rates)

#Population
   pop.times <<- pop.scenario$Year - 2014 #years
   L.N.y <<- pop.scenario$AdvancedRegion #millions
   L.S.y <<- pop.scenario$EmergingRegion #millions

#Calibrated parameters

   calib.params<<-c(
                 epsilon.N = round(2.269395134,4),  #round(2.288057866,4),
                 epsilon.S = round(3.166191929,4),  #round(3.012365201,4),
                 Gamma_re =  round(0.027114996,4), #round(0.029518117,4),
                 Gamma_ce=   round(0.006830185,4),  #round(0.007385199,4),
                 Eta_re.N =  round(0.761353468,4), #round(0.702797514,4),
                 Eta_ce.N =  round(0.942985835,4), #round(0.822026698,4),
                 Eta_re.S =  round(0.051594995,4), #round(0.043944150,4),
                 Eta_ce.S =  round(0.023885790,4), #round(0.022691853,4),
                 val.param.N = round(1.394925760,4), ##round(1.466980717,4),
                 val.param.S = round(1.499957895,4), #round(1.393887206,4),
                 lrng.re.N = 0.0,
                 lrng.ce.N = 0.0,
                 lrng.re.S = 0.0,
                 lrng.ce.S = 0.0,
                 alfa.N = round(0.316324021,4), #round(0.321840915,4),
                 alfa.S = round(0.299931522,4)  #round(0.321840915,4)
                 )
#Climate parameters

 climate.params<<-c(
                   qsi= climate.scenario$qsi,
                   Delta.S = climate.scenario$Delta.S,
                   Delta.Temp.Disaster = climate.scenario$Delta.Temp.Disaster,
                   Beta.Delta.Temp = climate.scenario$Beta.Delta.Temp,
                   CO2.base = climate.scenario$CO2.base,
                   CO2.Disaster = climate.scenario$CO2.Disaster,
                   DeltaTempConstant = climate.scenario$DeltaTempConstant,
                   S_0 = climate.scenario$S.0
                   )

#x<-c(0,0,0,0,0,0,0,0)
#Policy vector
#   policy.params<<-c(
#                  #Advanced Region
#                      ce.tax.N=round(x[1],3),
#                      Tec.subsidy.N = round(x[2],3),
#                      RD.subsidy.N = round(x[3],3),
#                      policy.half.life.N = round(x[4],4),
#                  #Emerging Region
#                      ce.tax.S=round(x[5],3),
#                      Tec.subsidy.S = round(x[6],3),
#                      RD.subsidy.S = round(x[7],3),
#                      policy.half.life.S = round(x[8],4)
#                    )
#
#Policy vector
   policy.params<<-c(
                  #Advanced Region
                      ce.tax.N = round(x[1],3),
                       Schedule.ce.tax.N = round(x[2],4),
                      Tec.subsidy.N = round(x[3],3),
                       Schedule.Tec.subsidy.N = round(x[4],4),
                      RD.subsidy.N = round(x[5],3),
                       Schedule.RD.subsidy.N = round(x[6],4),
                  #Emerging Region
                      ce.tax.S = round(x[7],3),
                       Schedule.ce.tax.S = round(x[8],4),
                      Tec.subsidy.S = round(x[9],3),
                       Schedule.Tec.subsidy.S = round(x[10],4),
                      RD.subsidy.S = round(x[11],3),
                       Schedule.RD.subsidy.S = round(x[12],4)
                    )

# Execute simulation
if (verbose==FALSE) {
                      Objective.Function.Value<-ediamMain(calib.params,verbose=FALSE)
                      return(Objective.Function.Value)
                    } else {
                      SimulData<-ediamMain(calib.params,verbose=TRUE)
                      return(SimulData)
                           }
}

#############################################################################################
# Use algorithm genoud to estimate parameters
#############################################################################################
#Optimization
#Load library for parallel computing
  library(snow)
#Specify numbers of cores available for calibration
  #nCore<-8
  nCore<-parallel::detectCores() - 1
#Define cluster
  cl <- makeSOCKcluster(names = rep('localhost',nCore))
  global.elements<-list("ediamPolicyEval","ediamMain","ediamEmpiricalParameters","ediamInitialConditions","EdiamEquations","Ediam","oil.scenario","pop.scenario","climate.scenario")
  clusterExport(cl,global.elements,envir=environment())

#Load libary genoud
 library(rgenoud)
#set seed for genetic optimization
 set.seed(55555)
#Execute the optimization
 genoud(ediamPolicyEval,max=TRUE,
       #nvars=8,
       nvars=12,
       starting.values = c( 0.5, #ce.tax.N
                             0.02, #Schedule.ce.tax.N
                            0.10, # Tec.subsidy.N
                             0.02, #Schedule.Tec.subsidy.N
                            1.0, # RD.Subsidy.N
                             0.02, #Schedule.RD.subsidy.N
                            #0.02, #policy.half.life.N
                            0.5, #ce.tax.S
                             0.02, #Schedule.ce.tax.S
                            0.10, # Tec.subsidy.S
                             0.02, #Schedule.Tec.subsidy.S
                            1.0, # RD.Subsidy.S
                             0.02 #Schedule.RD.subsidy.S
                            #0.02 #policy.half.life.S
                          ),
       pop.size=1000,
       Domains=matrix(c(
                        #inferior limits
                              0.0, #ce.tax.N
                               0.0, #Schedule.ce.tax.N
                              0.0, # Tec.subsidy.N
                               0.0, #Schedule.Tec.subsidy.N
                              0.0, # RD.Subsidy.N
                               0.0, #Schedule.RD.subsidy.N
                              #0.0, #policy.half.life.N
                              0.0, #ce.tax.S
                               0.0, #Schedule.ce.tax.S
                              0.0, # Tec.subsidy.S
                               0.0, #Schedule.Tec.subsidy.S
                              0.0, # RD.Subsidy.S
                               0.0, #Schedule.RD.subsidy.S
                              #0.0, #policy.half.life.S
                        #superior limits
                              1.0, #ce.tax.N
                               0.1, #Schedule.ce.tax.N
                              0.9, # Tec.subsidy.N
                               0.1, #Schedule.Tec.subsidy.N
                              4.0, # RD.Subsidy.N
                               0.1, #Schedule.RD.subsidy.N
                              #0.1, #policy.half.life.N
                              1.0, #ce.tax.S
                               0.1, #Schedule.ce.tax.S
                              0.9, # Tec.subsidy.S
                               0.1, #Schedule.Tec.subsidy.S
                              4.0, # RD.Subsidy.S
                               0.1 #Schedule.RD.subsidy.S
                              #0.1 #policy.half.life.S
                              ),
                            ncol=2),
       cluster=cl,
       print.level=1)

stopCluster(cl)

#Solution for climate scenario:"GFDL-ESM2G" , population scenario:"UN.Median.PI", timestep=1, lambda.S = 0.1443*4.0
#Value of objetive function: -41.2252
#solution
 x<-c(
                #Advanced Region
                 0.1759729267 ,  #ce.tax.N
                   0.0083342075,
                 0.8278753097,   #Tec.subsidy.N
                   0.0023607783,
                0.2428055360,   #RD.subsidy.N
                   0.0748504300,   #
                #Emerging Region
                 0.1915780205,   #ce.tax.S
                   -0.0007484287,
                 0.6971772584,  #Tec.subsidy.S
                   0.0091336108,
                 0.0002592633,   #RD.subsidy.S
                   0.0652135338
                  )

##Solution for climate scenario:"GFDL-ESM2G" , population scenario:"UN.Lower95.PI", timestep=5, lambda.S = 0.1443*4.0
#Value of objetive function: -10.92281
x<-c(0.260432451,0.019417471,0.833999113,0.002035506,0.540403727,0.080215080,0.180240385,0.004941565,0.696275779,0.001916819,0.001107793,0.010012071)

##Solution for climate scenario:"GFDL-ESM2G" , population scenario:"UN.Median.PI", timestep=5, lambda.S = 0.1443*4.0
#Value of objetive function: -10.9294
x<-c(0.1633620908,0.0043222848,0.8239402876,0.0015262409,0.6740002979,0.0749896069,0.1645699259,0.0018567055,0.7039641752,0.0022420558,0.0003846924,0.044151790)

# recupera la simulación completa con los parámetros óptimos
SimulData <- ediamPolicyEval(x, verbose = TRUE)

#cost of policy
 total.cost.N<-sum(SimulData$RealTecsubsidy.N*5)+sum(SimulData$RealRDsubsidy.N*5)
 total.cost.N
# Share on RD
   sum(SimulData$RealRDsubsidy.N*5)/total.cost.N
# Share on Technology
   sum(SimulData$RealTecsubsidy.N*5)/total.cost.N
#
 total.cost.S<-sum(SimulData$RealTecsubsidy.S*5)+sum(SimulData$RealRDsubsidy.S*5)
 total.cost.S
# Share on RD
   sum(SimulData$RealRDsubsidy.S*5)/total.cost.S
# Share on Technology
   sum(SimulData$RealTecsubsidy.S*5)/total.cost.S

#Share of GDP
  total.cost.N/sum(SimulData$GDP.N*5)
  total.cost.S/sum(SimulData$GDP.S*5)



x<-c(
               #Advanced Region
                0.0,  #ce.tax.N
                0.0,
                0.0,   #Tec.subsidy.N
                  0.0,
                0.0,   #RD.subsidy.N
                0.0,   #policy.half.life.N
               #Emerging Region
                0.0,   #ce.tax.S
                0.0,
                0.0,  #Tec.subsidy.S
                0.0,
                0.0,   #RD.subsidy.S
                0.0 #policy.half.life.S
                 )

#

ediamPolicyEval(x)

#first we would need to simulate the model

SimulData<-ediamPolicyEval(x,verbose=TRUE)

library(ggplot2)
library(reshape2)
library(ggpubr)

ediamPlot(SimulData)

#with new equations and parameters
#if you do nothing - > objective function: -44.46336;  if cost of damage is low: -43.14165

#if climate has no effect on you,  objective function; if cost of damage is low:  -42.68452

#if you keep temparature below 2 celsius: -41.84196 ; if cost of damage is low: -41.14351

#############################################################################################
# The following script uses the function ediamMSPE in combination with the genoud optimization
# algorithm to estimate the set of parameters that minimize MSPE with respect to historical data
# For robustness, the calibration attempt is executed in four different experimentes
# I.   Using historical data from 1995 to 2014 (i.e. monotonic increase in oil price) and no learning effects
# II.  Using historical data from 1971 to 2014 and no learning effects
# III. Using historical data from 1995 to 2014 (i.e. monotonic increase in oil price) and learning effects
# IV.  Using historical data from 1971 to 2014 and learning effects
#############################################################################################

#=============================================================================================
# Experiment I: Using historical data from 1995 to 2014 (i.e. monotonic increase in oil price) and no learning effects
#============================================================================================
#Load library for using HP filter
 library(mFilter)
#Read data file
  HistData<<-read.csv(file.path(dir.data, "HistData.csv"))
#Subset the data to the years that we have data for
#   HistData<<-subset(HistData,HistData$Time%in%c(1995:2014))
   HistData<<-subset(HistData,HistData$Time%in%c(1971:2014))
   HistData$Y_ce.h<-HistData[,"TPES"]*HistData[,"FOSSILTPES"] #MTOE
   HistData$Y_re.h<-HistData[,"TPES"]*HistData[,"RENTPES"] #MTOE
#Create NonOECD region
  HistDataNonOECD<-subset(HistData,HistData$Country=="World")
  HistDataNonOECD$Country<-"NonOECD"
  HistDataNonOECD$CountryID<-"NonOECD"
  HistDataNonOECD$TPES<-HistDataNonOECD$TPES-HistData$TPES[HistData$Country=="Memo: OECD Total"]
  HistDataNonOECD$Y_ce.h<-HistDataNonOECD$Y_ce.h-HistData$Y_ce.h[HistData$Country=="Memo: OECD Total"]
  HistDataNonOECD$Y_re.h<-HistDataNonOECD$Y_re.h-HistData$Y_re.h[HistData$Country=="Memo: OECD Total"]
  HistDataNonOECD$GDP<-HistDataNonOECD$GDP-HistData$GDP[HistData$Country=="Memo: OECD Total"]
  HistDataNonOECD$POP<-HistDataNonOECD$POP-HistData$POP[HistData$Country=="Memo: OECD Total"]
  HistDataNonOECD$OILTPES<-HistDataNonOECD$OILTPES-HistData$OILTPES[HistData$Country=="Memo: OECD Total"]
#Rbind
  HistData<-rbind(HistData,HistDataNonOECD)
#Additional operations
   HistData$c_oil<-HistData[,"OILTPES"]/(HistData[,"TPES"]*HistData[,"FOSSILTPES"])#[1]
   HistData$Roil.h<-HistData[,"OILTPES"] #MTOE
   HistData$Re.h<-HistData$Roil.h/HistData$c_oil #MTOE
   HistData$Price.Oil<- HistData$Price.Oil*7.1428571428571*1e6/1e9 # (billion USD per Mtoe) ; Asumming that 1 toe= 7.1428571428571 boe
   HistData$ReToGDP.h<-(HistData$Re.h*HistData$Price.Oil)/(HistData$GDP) # [1]
#hp filter to oil prices
   HistData$Price.Oil.hp<-rep(as.numeric(hpfilter(subset(HistData[,"Price.Oil"],HistData[,"Country"]=="World"),freq=100)$trend),length(unique(HistData$Country)))

#Indicate variables of interest for estimating MSPE
  ValVars<<-c("GDP.N","Y_ce.N","Y_re.N","GDP.S","Y_ce.S","Y_re.S")

#Load script for running the historical calibration
  source(file.path(dir.calib, "ediam_HistoricCalib.r"))
#Load script for determining initial conditions
  source(file.path(dir.calib, "ediam_InitialConditions.r"))
#Load EDIAM Equations
  source(file.path(dir.model, "ediam_Equations_discrete.r"))

# Specify ediamMSPE function for this experiment

ediamMSPE<-function(x,verbose=FALSE){
 calib.params<<-c(
                   epsilon.N = round(x[1],4),
                   epsilon.S = round(x[2],4),
                   Gamma_re =  round(x[3],4),
                   Gamma_ce=   round(x[4],4),
                   Eta_re.N =  round(x[5],4),
                   Eta_ce.N =  round(x[6],4),
                   Eta_re.S =  round(x[7],4),
                   Eta_ce.S =  round(x[8],4),
                   val.param.N = round(x[9],4),
                   val.param.S = round(x[10],4),
                   lrng.re.N = 0.0, #Note: no experience accumulation effects considered in this experiment
                   lrng.ce.N = 0.0, #Note: no experience accumulation effects considered in this experiment
                   lrng.re.S = 0.0, #Note: no experience accumulation effects considered in this experiment
                   lrng.ce.S = 0.0, #Note: no experience accumulation effects considered in this experiment
                   alfa.N = round(x[11],4),
                   alfa.S = round(x[12],4)
                   )

#historic data with hp
  Y_re.Nh <<-as.numeric(hpfilter(subset(HistData[,"Y_re.h"],HistData[,"Country"]=="Memo: OECD Total"),freq=100)$trend) #Mtoe
  Y_ce.Nh <<-as.numeric(hpfilter(subset(HistData[,"Y_ce.h"],HistData[,"Country"]=="Memo: OECD Total"),freq=100)$trend) #Mtoe
  Y_re.Sh <<-as.numeric(hpfilter(subset(HistData[,"Y_re.h"],HistData[,"Country"]=="NonOECD"),freq=100)$trend) #Mtoe
  Y_ce.Sh <<-as.numeric(hpfilter(subset(HistData[,"Y_ce.h"],HistData[,"Country"]=="NonOECD"),freq=100)$trend) #Mtoe

#Oil Supply
  Re.Nh<<-subset(HistData[,"Re.h"],HistData[,"Country"]=="Memo: OECD Total") # Mtoe
  Re.Sh<<-subset(HistData[,"Re.h"],HistData[,"Country"]=="NonOECD") #Mtoe
  ReToGDP.Nh<<-subset(HistData[,"ReToGDP.h"],HistData[,"Country"]=="Memo: OECD Total") # [1]
  ReToGDP.Sh<<-subset(HistData[,"ReToGDP.h"],HistData[,"Country"]=="NonOECD")

#Oil prices
   Price.oil.y<<- subset(HistData[,"Price.Oil.hp"],HistData[,"Country"]=="World")

#GDP with hp
    GDP.Nh<<-as.numeric(hpfilter(subset(HistData[,"GDP"],HistData[,"Country"]=="Memo: OECD Total"),freq=100)$trend) #(billion 2010 USD using exchange rates)
    GDP.Sh<<-as.numeric(hpfilter(subset(HistData[,"GDP"],HistData[,"Country"]=="NonOECD"),freq=100)$trend) #(billion 2010 USD using exchange rates)

#Population
   L.N.y <<-subset(HistData[,"POP"],HistData[,"Country"]=="Memo: OECD Total") #millions
   L.S.y <<-subset(HistData[,"POP"],HistData[,"Country"]=="NonOECD") #millions

#Simulate model
   SimulData<-ediamCalib(calib.params,verbose=TRUE)
#Estimate mean square differences
  # MSPEp1<-((SimulData[,ValVars[1]]-GDP.Nh)/GDP.Nh)^2
  # MSPEp1.1<-((SimulData[,ValVars[2]]-Y_ce.Nh)/Y_ce.Nh)^2
  # MSPEp1.2<-((SimulData[,ValVars[3]]-Y_re.Nh)/Y_re.Nh)^2
  # MSPEp2<-((SimulData[,ValVars[4]]-GDP.Sh)/GDP.Sh)^2
  # MSPEp2.1<-((SimulData[,ValVars[5]]-Y_ce.Sh)/Y_ce.Sh)^2
  # MSPEp2.2<-((SimulData[,ValVars[6]]-Y_re.Sh)/Y_re.Sh)^2
   MSPEp1<-abs(SimulData[,ValVars[1]]-GDP.Nh)/GDP.Nh
   MSPEp1.1<-abs(SimulData[,ValVars[2]]-Y_ce.Nh)/Y_ce.Nh
   MSPEp1.2<-abs(SimulData[,ValVars[3]]-Y_re.Nh)/Y_re.Nh
   MSPEp2<-abs(SimulData[,ValVars[4]]-GDP.Sh)/GDP.Sh
   MSPEp2.1<-abs(SimulData[,ValVars[5]]-Y_ce.Sh)/Y_ce.Sh
   MSPEp2.2<-abs(SimulData[,ValVars[6]]-Y_re.Sh)/Y_re.Sh
#estimate mean MSE
   MSPE<-mean(c(mean(MSPEp1),mean(MSPEp1.1),mean(MSPEp1.2),mean(MSPEp2),mean(MSPEp2.1),mean(MSPEp2.2)))
#check constraints
#   tolerance<-0.011
   tolerance<-0.11
   MSPE<- ifelse(mean(MSPEp1)<=tolerance,
               ifelse(mean(MSPEp1.1)<=tolerance,
                       ifelse(mean(MSPEp1.2)<=tolerance,
                                ifelse(mean(MSPEp2)<=tolerance,
                                     ifelse(mean(MSPEp2.1)<=tolerance,
                                          ifelse(mean(MSPEp2.2)<=tolerance,MSPE,1e3)
                                                                 ,1e3)
                                                        ,1e3)
                                                    ,1e3)
                                            ,1e3)
                                    ,1e3)
# Concatenate results
  out<-data.frame(v1 = mean(MSPEp1),
                  v2 = mean(MSPEp1.1),
                  v3 = mean(MSPEp1.2),
                  v4 = mean(MSPEp2),
                  v5 = mean(MSPEp2.1),
                  v6 = mean(MSPEp2.2),
                  v7 = MSPE )

  colnames(out)<-paste("MSPE",c(ValVars[1],ValVars[2],ValVars[3],ValVars[4],ValVars[5],ValVars[6],"all"),sep=".")
if (verbose == TRUE)
{ return(out)
} else {
        return(as.numeric(out$MSPE.all))
       }

}

#############################################################################################
# Use algorithm genoud to estimate parameters
#############################################################################################
#Optimization
#Load library for parallel computing
  library(snow)
#Specify numbers of cores available for calibration
  nCore<-parallel::detectCores() - 1
#Define cluster
  cl <- makeSOCKcluster(names = rep('localhost',nCore))
  global.elements<-list("ediamMSPE","ediamCalib","ediamEmpiricalParameters","ediamInitialConditions","EdiamEquations","Ediam","HistData","ValVars","dir.model","hpfilter")
  clusterExport(cl,global.elements,envir=environment())

#Load libary genoud
 library(rgenoud)
#set seed for genetic optimization
 set.seed(55555)
#Execute the optimization
 genoud(ediamMSPE,max=FALSE,
       nvars=12,
       starting.values = c( 2.27,  #epsilon.N
                            2.91, #epsilon.S
                            0.024, #Gamma_re
                            0.006, #Gamma_ce
                            0.484, #Eta_re.N
                            0.64, #Eta_ce.N
                            0.052, #Eta_re.S
                            0.026, #, #Eta_ce.S
                            1.0,   #val.param.N
                            1.0,   #val.param.S
                          #  0.000 , #lrng.re.N
                          #  0.000 ,#, #lrng.ce.N
                          #  0.000 , #lrng.re.N
                          #  0.000 , #lrng.ce.S
                            0.333, #alfa.N
                            0.333 #alfa.S
                          ),
       pop.size=10000,
       Domains=matrix(c(#inferior limits
                            1.5,  #epsilon.N
                            1.5, #epsilon.S
                            0.001, #Gamma_re
                            0.001, #Gamma_ce
                            0.001, #Eta_re.N
                            0.001, #Eta_ce.N
                            0.001, #Eta_re.S
                            0.001, #Eta_ce.S
                            0.85, #val.param.N
                            0.85, #val.param.S
                        #    0.000 , #lrng.re.N
                        #    0.000 , #lrng.ce.N
                        #    0.000 , #lrng.re.S
                        #    0.000 , #lrng.ce.S,
                            0.15, #alfa.N
                            0.15, #alfa.S
                        #superior limits
                            10,  #epsilon.N
                            10, #epsilon.S
                            0.3, #Gamma_re
                            0.3, #Gamma_ce
                            0.95, #Eta_re.N
                            0.95, #Eta_ce.N
                            0.95, #Eta_re.S
                            0.95, #Eta_ce.S
                            1.5 , #val.param.N
                            1.5 , #val.param.S
                          #  1.0 , #lrng.re.N
                          #  1.0 ,#, #lrng.ce.N
                          #  1.0 , #lrng.re.S
                          #  1.0 , #lrng.ce.S
                            0.5, #alfa.N
                            0.5 #alfa.S
                               ),
                            ncol=2),
       cluster=cl,
       print.level=1)

stopCluster(cl)

#The solution for this experiment is the following
#MSPE<-0.004254538
x<-c(2.27000000,
     2.92492559,
     0.02571846,
     0.00600000,
     0.66104448,
     0.88082188,
     0.05200000,
     0.02600000,
     1.20210188,
     1.49992306,
     0.31991824,
     0.30724045
   )

#MAPE<-0.04265314 prefered solution *
x<-c( 2.288057866,
      3.012365201,
      0.029518117,
      0.007385199,
      0.702797514,
      0.822026698,
      0.043944150,
      0.022691853,
      1.466980717,
      1.393887206,
      0.321840915,
      0.322528624
    )

# verify solution
 ediamMSPE(x, verbose = TRUE)
#save calibrated parameters
  calib.params<<-c(
                  epsilon.N = round(x[1],4),
                  epsilon.S = round(x[2],4),
                  Gamma_re =  round(x[3],4),
                  Gamma_ce=   round(x[4],4),
                  Eta_re.N =  round(x[5],4),
                  Eta_ce.N =  round(x[6],4),
                  Eta_re.S =  round(x[7],4),
                  Eta_ce.S =  round(x[8],4),
                  val.param.N = round(x[9],4),
                  val.param.S = round(x[10],4),
                  lrng.re.N = 0.0, #Note: no experience accumulation effects considered in this experiment
                  lrng.ce.N = 0.0, #Note: no experience accumulation effects considered in this experiment
                  lrng.re.S = 0.0, #Note: no experience accumulation effects considered in this experiment
                  lrng.ce.S = 0.0, #Note: no experience accumulation effects considered in this experiment
                  alfa.N = round(x[11],4),
                  alfa.S = round(x[12],4)
                  )

#Compare graphically solution found vs historical record using the following routine
  Y_re.Nh <-as.numeric(hpfilter(subset(HistData[,"Y_re.h"],HistData[,"Country"]=="Memo: OECD Total"),freq=100)$trend) #Mtoe
  Y_ce.Nh <-as.numeric(hpfilter(subset(HistData[,"Y_ce.h"],HistData[,"Country"]=="Memo: OECD Total"),freq=100)$trend) #Mtoe
  Y_re.Sh <-as.numeric(hpfilter(subset(HistData[,"Y_re.h"],HistData[,"Country"]=="NonOECD"),freq=100)$trend) #Mtoe
  Y_ce.Sh <-as.numeric(hpfilter(subset(HistData[,"Y_ce.h"],HistData[,"Country"]=="NonOECD"),freq=100)$trend) #Mtoe
  Re.Nh<-subset(HistData[,"Re.h"],HistData[,"Country"]=="Memo: OECD Total") # Mtoe
  Re.Sh<-subset(HistData[,"Re.h"],HistData[,"Country"]=="NonOECD") #Mtoe
  ReToGDP.Nh<-subset(HistData[,"ReToGDP.h"],HistData[,"Country"]=="Memo: OECD Total") # [1]
  ReToGDP.Sh<-subset(HistData[,"ReToGDP.h"],HistData[,"Country"]=="NonOECD")
  Price.oil.y<- subset(HistData[,"Price.Oil.hp"],HistData[,"Country"]=="World")
  GDP.Nh<-as.numeric(hpfilter(subset(HistData[,"GDP"],HistData[,"Country"]=="Memo: OECD Total"),freq=100)$trend) #(billion 2010 USD using exchange rates)
  GDP.Sh<-as.numeric(hpfilter(subset(HistData[,"GDP"],HistData[,"Country"]=="NonOECD"),freq=100)$trend) #(billion 2010 USD using exchange rates)
  L.N.y <-subset(HistData[,"POP"],HistData[,"Country"]=="Memo: OECD Total") #millions
  L.S.y <-subset(HistData[,"POP"],HistData[,"Country"]=="NonOECD") #millions

#Simulate with calibrated parameters
  SimulData<-ediamCalib(calib.params,verbose=TRUE)

#Reshape simulated data into usable format for comparison
 times<-seq(min(HistData$Time),max(HistData$Time))
 SimulData<-data.frame( Time = rep(times,7),
           Value = c(SimulData[,"GDP.N"],SimulData[,"GDP.S"], SimulData[,"Y_re.N"], SimulData[,"Y_ce.N"],SimulData[,"Y_re.S"],SimulData[,"Y_ce.S"],SimulData[,"Price.oil"]),
           Variable= rep(c("GDP.OECD","GDP.NonOECD","Y_re.OECD","Y_ce.OECD","Y_re.NonOECD","Y_ce.NonOECD","Price.Oil"),each=length(times)),
           Data.type="Simulation")

#Reshape historcial data into usable format for comparison
 HistoricData.hp<-data.frame( Time = rep(times,7),
           Value = c(GDP.Nh,GDP.Sh,Y_re.Nh,Y_ce.Nh,Y_re.Sh,Y_ce.Sh,Price.oil.y),
           Variable= rep(c("GDP.OECD","GDP.NonOECD","Y_re.OECD","Y_ce.OECD","Y_re.NonOECD","Y_ce.NonOECD","Price.Oil"),each=length(times)),
           Data.type="Historic")

CalibFitness<-rbind(SimulData,HistoricData.hp)

library(ggplot2)
#GDP
calib.plot <- ggplot(data=subset(CalibFitness,CalibFitness$Variable%in%c("GDP.OECD","GDP.NonOECD")), aes(x=Time, y=Value, linetype = Data.type))
calib.plot + geom_line()+scale_linetype_manual(values=c(2,1)) + facet_wrap(~ Variable)

#FOSSILTPES
calib.plot <- ggplot(data=subset(CalibFitness,CalibFitness$Variable%in%c("Y_ce.OECD","Y_ce.NonOECD")), aes(x=Time, y=Value, linetype = Data.type))
calib.plot + geom_line()+scale_linetype_manual(values=c(2,1)) + facet_wrap(~ Variable)

#RENTPES
calib.plot <- ggplot(data=subset(CalibFitness,CalibFitness$Variable%in%c("Y_re.OECD","Y_re.NonOECD")), aes(x=Time, y=Value, linetype = Data.type))
calib.plot + geom_line()+scale_linetype_manual(values=c(2,1)) + facet_wrap(~ Variable)



#=============================================================================================
# Experiment II: Using historical data from 1971 to 2014 (i.e. monotonic increase in oil price) and no learning effects
#============================================================================================



#MAPE<-0.08010907
x<-c(1.123808,
     1.359642,
     0.01500615,
     0.01047485,
     0.40957,
     0.5420697,
     0.01260745,
     0.01278657,
     1.082811,
     1.063265,
     0.07523035,
     0.1943299
    )

# verify solution
 ediamMSPE(x, verbose = TRUE)
#save calibrated parameters
  calib.params<<-c(
                  epsilon.N = round(x[1],4),
                  epsilon.S = round(x[2],4),
                  Gamma_re =  round(x[3],4),
                  Gamma_ce=   round(x[4],4),
                  Eta_re.N =  round(x[5],4),
                  Eta_ce.N =  round(x[6],4),
                  Eta_re.S =  round(x[7],4),
                  Eta_ce.S =  round(x[8],4),
                  val.param.N = round(x[9],4),
                  val.param.S = round(x[10],4),
                  lrng.re.N = 0.0, #Note: no experience accumulation effects considered in this experiment
                  lrng.ce.N = 0.0, #Note: no experience accumulation effects considered in this experiment
                  lrng.re.S = 0.0, #Note: no experience accumulation effects considered in this experiment
                  lrng.ce.S = 0.0, #Note: no experience accumulation effects considered in this experiment
                  alfa.N = round(x[11],4),
                  alfa.S = round(x[12],4)
                  )

#Compare graphically solution found vs historical record using the following routine
  Y_re.Nh <-as.numeric(hpfilter(subset(HistData[,"Y_re.h"],HistData[,"Country"]=="Memo: OECD Total"),freq=100)$trend) #Mtoe
  Y_ce.Nh <-as.numeric(hpfilter(subset(HistData[,"Y_ce.h"],HistData[,"Country"]=="Memo: OECD Total"),freq=100)$trend) #Mtoe
  Y_re.Sh <-as.numeric(hpfilter(subset(HistData[,"Y_re.h"],HistData[,"Country"]=="NonOECD"),freq=100)$trend) #Mtoe
  Y_ce.Sh <-as.numeric(hpfilter(subset(HistData[,"Y_ce.h"],HistData[,"Country"]=="NonOECD"),freq=100)$trend) #Mtoe
  Re.Nh<-subset(HistData[,"Re.h"],HistData[,"Country"]=="Memo: OECD Total") # Mtoe
  Re.Sh<-subset(HistData[,"Re.h"],HistData[,"Country"]=="NonOECD") #Mtoe
  ReToGDP.Nh<-subset(HistData[,"ReToGDP.h"],HistData[,"Country"]=="Memo: OECD Total") # [1]
  ReToGDP.Sh<-subset(HistData[,"ReToGDP.h"],HistData[,"Country"]=="NonOECD")
  Price.oil.y<- subset(HistData[,"Price.Oil.hp"],HistData[,"Country"]=="World")
  GDP.Nh<-as.numeric(hpfilter(subset(HistData[,"GDP"],HistData[,"Country"]=="Memo: OECD Total"),freq=100)$trend) #(billion 2010 USD using exchange rates)
  GDP.Sh<-as.numeric(hpfilter(subset(HistData[,"GDP"],HistData[,"Country"]=="NonOECD"),freq=100)$trend) #(billion 2010 USD using exchange rates)
  L.N.y <-subset(HistData[,"POP"],HistData[,"Country"]=="Memo: OECD Total") #millions
  L.S.y <-subset(HistData[,"POP"],HistData[,"Country"]=="NonOECD") #millions

#Simulate with calibrated parameters
  SimulData<-ediamCalib(calib.params,verbose=TRUE)

#Reshape simulated data into usable format for comparison
 times<-seq(min(HistData$Time),max(HistData$Time))
 SimulData<-data.frame( Time = rep(times,7),
           Value = c(SimulData[,"GDP.N"],SimulData[,"GDP.S"], SimulData[,"Y_re.N"], SimulData[,"Y_ce.N"],SimulData[,"Y_re.S"],SimulData[,"Y_ce.S"],SimulData[,"Price.oil"]),
           Variable= rep(c("GDP.OECD","GDP.NonOECD","Y_re.OECD","Y_ce.OECD","Y_re.NonOECD","Y_ce.NonOECD","Price.Oil"),each=length(times)),
           Data.type="Simulation")

#Reshape historcial data into usable format for comparison
 HistoricData.hp<-data.frame( Time = rep(times,7),
           Value = c(GDP.Nh,GDP.Sh,Y_re.Nh,Y_ce.Nh,Y_re.Sh,Y_ce.Sh,Price.oil.y),
           Variable= rep(c("GDP.OECD","GDP.NonOECD","Y_re.OECD","Y_ce.OECD","Y_re.NonOECD","Y_ce.NonOECD","Price.Oil"),each=length(times)),
           Data.type="Historic")

CalibFitness<-rbind(SimulData,HistoricData.hp)

library(ggplot2)
#GDP
calib.plot <- ggplot(data=subset(CalibFitness,CalibFitness$Variable%in%c("GDP.OECD","GDP.NonOECD")), aes(x=Time, y=Value, linetype = Data.type))
calib.plot + geom_line()+scale_linetype_manual(values=c(2,1)) + facet_wrap(~ Variable)

#FOSSILTPES
calib.plot <- ggplot(data=subset(CalibFitness,CalibFitness$Variable%in%c("Y_ce.OECD","Y_ce.NonOECD")), aes(x=Time, y=Value, linetype = Data.type))
calib.plot + geom_line()+scale_linetype_manual(values=c(2,1)) + facet_wrap(~ Variable)

#RENTPES
calib.plot <- ggplot(data=subset(CalibFitness,CalibFitness$Variable%in%c("Y_re.OECD","Y_re.NonOECD")), aes(x=Time, y=Value, linetype = Data.type))
calib.plot + geom_line()+scale_linetype_manual(values=c(2,1)) + facet_wrap(~ Variable)



#=================
#alfa parameters, all data
#with 10,000 pop
#MSE<-0.01488052

x<-c(1.627149819,
     1.814917316,
     0.012922016,
     0.004417407,
     0.310772866,
     0.660008525,
     0.031857169,
     0.019505122,
     1.145331810,
     1.024535191,
     0.327010721,
     0.205050104
   )

#alfa parameters, subset of data, learning Effects
#MSE<-0.004193859

x<-c(2.27000000,
     2.91000000,
     0.02453531,
     0.00600000,
     0.52577957,
     0.69132343,
     0.05418693,
     0.02600000,
     1.00000000,
     1.49895532,
     0.00000000,
     0.00000000,
     0.10126171,
     0.27304493,
     0.33193600,
     0.31014238
   )

#alfa parameters, subset of data, learning Effects
#MSE<-0.004193859
#tolerance<-0.011

x<-c(2.270008e+00,
     2.879429e+00,
     2.490752e-02,
     5.687152e-03,
     4.840018e-01,
     7.773862e-01,
     5.260620e-02 ,
     2.531791e-02 ,
     8.806384e-01 ,
     1.425176e+00 ,
     1.028267e-01 ,
     -8.921408e-07,
     1.384415e-01,
     6.128041e-01,
     3.371347e-01,
     3.223382e-01
   )
##alfa parameters, all data, learning Effects
#MSE<-0.01234146

x<-c(1.537581177,
     1.840202353,
     0.005459022,
     0.003417549,
     0.176664977,
     0.559690277,
     0.003501534,
     0.004122687,
     1.075556602,
     1.193103218,
     0.472986492,
     0.041228237,
     0.310514641,
     0.998759662,
     0.408460563,
     0.386325317
   )




source( file.path(dir.calib, "ediam_HistoricCalib.r") )
#first show MSE
ediamMSPE(x, verbose = TRUE)
#Compare visually the fittness of the solution
#
calib.params<<-c(
                  epsilon.N = round(x[1],3),
                  epsilon.S = round(x[2],3),
                  Gamma_re =  round(x[3],3),
                  Gamma_ce=   round(x[4],3),
                  Eta_re.N =  round(x[5],3),
                  Eta_ce.N =  round(x[6],3),
                  Eta_re.S =  round(x[7],3),
                  Eta_ce.S =  round(x[8],3),
                  val.param.N = round(x[9],3),
                  val.param.S = round(x[10],3),
                  pi.N = 1.0, #round(x[10],3),
                  pi.S = 1.0, #round(x[11],3),
                  lrng.re.N = round(x[11],3),
                  lrng.ce.N = round(x[12],3), #,
                  lrng.re.S = round(x[13],3),
                  lrng.ce.S = round(x[14],3),
                  cf.re.N = 1.0, #round(x[9],3),
                  cf.ce.N = 1.0, #round(x[10],3),
                  cf.re.S = 1.0, #round(x[11],3),
                  cf.ce.S = 1.0, #round(x[12],3)
                  alfa.N = round(x[15],3),
                  alfa.S = round(x[16],3)
                  )

 SimulData<-ediamCalib(calib.params, verbose=TRUE)
 SimulData$GDPpc.OECD<-SimulData$GDP.N/SimulData$L.N
 SimulData$GDPpc.NonOECD<-SimulData$GDP.S/SimulData$L.S


#design the effects of experience
#  v<-SimulData$Xtech_ce.N
#  Experience <- sapply(1:(length(v)),function(x){sum(v[1:x])})
#  lrng.re<-1.0 #range from 0 to 1.0
#  cost <- exp(-1* lrng.re* 1e-2 *(Experience-v[1])/v[1])

#Reshape simulated data into usable format for comparison
 times<-seq(min(HistData$Time),max(HistData$Time))
 SimLong<-data.frame( Time = rep(times,11),
            Value = c(SimulData[,"GDP.N"],SimulData[,"GDP.S"], SimulData[,"Y_re.N"], SimulData[,"Y_ce.N"],SimulData[,"Y_re.S"],SimulData[,"Y_ce.S"],
                      SimulData[,"L.N"],SimulData[,"L.S"],SimulData[,"Price.oil"],SimulData[,"GDPpc.OECD"],SimulData[,"GDPpc.NonOECD"]),
            Variable= rep(c("GDP.OECD","GDP.NonOECD","Y_re.OECD","Y_ce.OECD","Y_re.NonOECD","Y_ce.NonOECD",
                            "POP.OECD","POP.NonOECD","Price.Oil","GDPpc.OECD","GDPpc.NonOECD"),each=length(times)),
            Data.type="Simulation")

#Reshape historcial data into usable format for comparison
 HistoricData.hp<-data.frame( Time = rep(times,11),
            Value = c(hpfilter(HistData[HistData$Country=="Memo: OECD Total","GDP"],freq=100)$trend,
                      hpfilter(HistData[HistData$Country=="NonOECD","GDP"],freq=100)$trend,
                      hpfilter(HistData[HistData$Country=="Memo: OECD Total","Y_re.h"],freq=100)$trend,
                      hpfilter(HistData[HistData$Country=="Memo: OECD Total","Y_ce.h"],freq=100)$trend,
                      hpfilter(HistData[HistData$Country=="NonOECD","Y_re.h"],freq=100)$trend,
                      hpfilter(HistData[HistData$Country=="NonOECD","Y_ce.h"],freq=100)$trend,
                      hpfilter(HistData[HistData$Country=="Memo: OECD Total","POP"],freq=100)$trend,
                      hpfilter(HistData[HistData$Country=="NonOECD","POP"],freq=100)$trend,
                      hpfilter(HistData[HistData$Country=="Memo: OECD Total","Price.Oil"],freq=100)$trend,
                      hpfilter(HistData[HistData$Country=="Memo: OECD Total","GDP"]/HistData[HistData$Country=="Memo: OECD Total","POP"],freq=100)$trend,
                      hpfilter(HistData[HistData$Country=="NonOECD","GDP"]/HistData[HistData$Country=="NonOECD","POP"],freq=100)$trend),
            Variable= rep(c("GDP.OECD","GDP.NonOECD","Y_re.OECD","Y_ce.OECD","Y_re.NonOECD","Y_ce.NonOECD",
                           "POP.OECD","POP.NonOECD","Price.Oil","GDPpc.OECD","GDPpc.NonOECD"),each=length(times)),
            Data.type="Historic")

#use hp on all historic data
HistoricData<-data.frame( Time = rep(times,11),
           Value = c(HistData[HistData$Country=="Memo: OECD Total","GDP"],
                     HistData[HistData$Country=="NonOECD","GDP"],
                     HistData[HistData$Country=="Memo: OECD Total","Y_re.h"],
                     HistData[HistData$Country=="Memo: OECD Total","Y_ce.h"],
                     HistData[HistData$Country=="NonOECD","Y_re.h"],
                     HistData[HistData$Country=="NonOECD","Y_ce.h"],
                     HistData[HistData$Country=="Memo: OECD Total","POP"],
                     HistData[HistData$Country=="NonOECD","POP"],
                     HistData[HistData$Country=="Memo: OECD Total","Price.Oil"],
                     HistData[HistData$Country=="Memo: OECD Total","GDP"]/HistData[HistData$Country=="Memo: OECD Total","POP"],
                     HistData[HistData$Country=="NonOECD","GDP"]/HistData[HistData$Country=="NonOECD","POP"]),
           Variable= rep(c("GDP.OECD","GDP.NonOECD","Y_re.OECD","Y_ce.OECD","Y_re.NonOECD","Y_ce.NonOECD",
                          "POP.OECD","POP.NonOECD","Price.Oil","GDPpc.OECD","GDPpc.NonOECD"),each=length(times)),
           Data.type="Historic")

#Calibration fitness
# CalibFitness<-rbind(SimulData,HistoricData)
 CalibFitness<-rbind(SimLong,HistoricData.hp)
library(ggplot2)
#GDP
 calib.plot <- ggplot(data=subset(CalibFitness,CalibFitness$Variable%in%c("GDP.OECD","GDP.NonOECD")), aes(x=Time, y=Value, linetype = Data.type))
 calib.plot + geom_line()+scale_linetype_manual(values=c(2,1)) + facet_wrap(~ Variable)

#GDP per capita
  calib.plot <- ggplot(data=subset(CalibFitness,CalibFitness$Variable%in%c("GDPpc.OECD","GDPpc.NonOECD")), aes(x=Time, y=Value, linetype = Data.type))
  calib.plot + geom_line()+scale_linetype_manual(values=c(2,1)) + facet_wrap(~ Variable)

#FOSSILTPES
calib.plot <- ggplot(data=subset(CalibFitness,CalibFitness$Variable%in%c("Y_ce.OECD","Y_ce.NonOECD")), aes(x=Time, y=Value, linetype = Data.type))
calib.plot + geom_line()+scale_linetype_manual(values=c(2,1)) + facet_wrap(~ Variable)

#RENTPES
calib.plot <- ggplot(data=subset(CalibFitness,CalibFitness$Variable%in%c("Y_re.OECD","Y_re.NonOECD")), aes(x=Time, y=Value, linetype = Data.type))
calib.plot + geom_line()+scale_linetype_manual(values=c(2,1)) + facet_wrap(~ Variable)


#notes for tomorrow
  # 1) check why I cannot subset the time series
  # 2) check how the errors look if you see GDP per capita, energy primary supply in rates-> errors look the same
  # 3) check how the oil and population plots look like -it may be necessary to filter the oil data series-

#next steps
  # try calibration with the utility parameter
  # try calibration with learning effects
  # try calibration with tolerance value

#then,it is done


# Population
 calib.plot <- ggplot(data=subset(CalibFitness,CalibFitness$Variable%in%c("POP.OECD","POP.NonOECD")), aes(x=Time, y=Value, linetype = Data.type))
 calib.plot + geom_line()+scale_linetype_manual(values=c(2,1)) + facet_wrap(~ Variable)
#oil price
 calib.plot <- ggplot(data=subset(CalibFitness,CalibFitness$Variable=="Price.Oil"), aes(x=Time, y=Value, linetype = Data.type))
 calib.plot + geom_line()+scale_linetype_manual(values=c(2,1)) + facet_wrap(~ Variable)
# hp oil price
 library(mFilter)
 Price.oil.hp<- hpfilter(subset(HistData[,"Price.Oil"],HistData[,"Country"]=="World"),freq=50)
 Price.oil.hp<-data.frame(Time=c(1971:2014),Price.oil.hp=as.numeric(Price.oil.hp$trend))
 calib.plot <- ggplot(data=Price.oil.hp, aes(x=Time, y=Price.oil.hp))
 calib.plot + geom_line()


#Show transition plot

#OECD
 TPES.OECD<-subset(SimLong,SimLong$Variable%in%c("Y_re.OECD","Y_ce.OECD"))
 TPES.OECD<-aggregate(Value ~ Time, data = TPES.OECD, FUN = sum)
 TPES.OECD$Region<-'OECD'

#NonOECD
 TPES.NonOECD<-subset(SimLong,SimLong$Variable%in%c("Y_re.NonOECD","Y_ce.NonOECD"))
 TPES.NonOECD<-aggregate(Value ~ Time, data = TPES.NonOECD, FUN = sum)
 TPES.NonOECD$Region<-'NonOECD'

#rbind
 Data<-rbind(TPES.OECD,TPES.NonOECD)
 Data$Type<-'Simulation'

#Prepare Historic Data
 TPES.Hist<-subset(HistData,HistData$Country%in%c("Memo: OECD Total","NonOECD"))
 TPES.Hist$Value<-(TPES.Hist$FOSSILTPES+TPES.Hist$RENTPES)*TPES.Hist$TPES
 TPES.Hist$Region<-ifelse(TPES.Hist$Country=='Memo: OECD Total','OECD','NonOECD')
 TPES.Hist$Type<-'Historic'
 TPES.Hist<-TPES.Hist[,c('Time','Region','Value','Type')]

#rbind all
 Data<-rbind(Data,TPES.Hist)

calib.plot <- ggplot(data=Data, aes(x=Time, y=Value, colour = Region))
calib.plot + geom_line() + facet_wrap(~ Type)


#designing the value function of entrepreneurs
 SimulData$Profits_re.N
 SimulData$Profits_ce.N

 x<-c(SimulData$Profits_re.N,SimulData$Profits_ce.N)
 x<-c(1:1000)
#type 1

 Value.t1<-data.frame(x=x)
 #Value.t1$Util<-log(0.2*Value.t1$x)
 Value.t1$Util<-log(Value.t1$x^1.1)
 Value.t1$type<-"type1"
#
#type 2
 Value.t2<-data.frame(x=x)
 #Value.t2$Util<-log(0.4*Value.t2$x)
 Value.t2$Util<-log(Value.t2$x^1.2)
 Value.t2$type<-"type2"
#
#type 3
 Value.t3<-data.frame(x=x)
 #Value.t3$Util<-log(0.6*Value.t3$x)
 Value.t3$Util<-log(Value.t3$x^1.3)
 Value.t3$type<-"type3"
#
#type 4
 Value.t4<-data.frame(x=x)
 #Value.t4$Util<-log(0.8*Value.t4$x)
 Value.t4$Util<-log(Value.t4$x^1.4)
 Value.t4$type<-"type4"
#
#type 5
 Value.t5<-data.frame(x=x)
 #Value.t5$Util<-log(1.0*Value.t5$x)
 Value.t5$Util<-log(Value.t5$x^1.5)
 Value.t5$type<-"type5"
#
#type 6
 Value.t6<-data.frame(x=x)
 #Value.t6$Util<-log(1.2*Value.t6$x)
 Value.t6$Util<-log(Value.t6$x^1.6)
 Value.t6$type<-"type6"
#
#type 7
 Value.t7<-data.frame(x=x)
 #Value.t7$Util<-log(1.4*Value.t7$x)
 Value.t7$Util<-log(Value.t7$x^1.7)
 Value.t7$type<-"type7"
#
#type 8
 Value.t8<-data.frame(x=x)
 #Value.t8$Util<-log(1.6*Value.t8$x)
 Value.t8$Util<-log(Value.t8$x^1.8)
 Value.t8$type<-"type8"
#
#type 9
 Value.t9<-data.frame(x=x)
 #Value.t9$Util<-log(1.8*Value.t9$x)
 Value.t9$Util<-log(Value.t9$x^1.9)
 Value.t9$type<-"type9"
#
#type 10
 Value.t10<-data.frame(x=x)
# Value.t10$Util<-log(2.0*Value.t10$x)
 Value.t10$Util<-log(Value.t10$x^1.10)
 Value.t10$type<-"type10"

#

Value.S<-rbind(Value.t1,Value.t2,Value.t3,Value.t4,Value.t5,Value.t6,Value.t7,Value.t8,Value.t9,Value.t10)

ggplot(data=Value.S, aes(x=x, y=Util,linetype=type))+geom_line()

#now let´s simulate this

  test<-data.frame( RE=SimulData$Profits_re.N,
                    CE=SimulData$Profits_ce.N,
                    Value.RE=log((SimulData$Profits_re.N/SimulData$Profits_re.N[1])^2.5),
                    Value.CE=log((SimulData$Profits_ce.N/SimulData$Profits_re.N[1])^2.5)
                  )
 test$s.re<- exp(test$Value.RE)/(exp(test$Value.CE)+exp(test$Value.CE))



#error analysis

library(lhs)
set.seed(5000)
# sample.size<-100000
sample.size<-1000
params.names<-c("epsilon.N",
                "epsilon.S",
                "Gamma_re",
                "Gamma_ce",
                "Eta_re.N",
                "Eta_ce.N",
                "Eta_re.S",
                "Eta_ce.S",
                "val.param.N",
                "val.param.S",
                "alfa.N",
                "alfa.S")
#
lhs.sample<-data.frame(randomLHS(sample.size,length(params.names)),Run.ID=1:sample.size)
colnames(lhs.sample)<-c(params.names,"Run.ID")
#Define range of values
lhs.sample$epsilon.N<-qunif(lhs.sample$epsilon.N,1.5,10)
lhs.sample$epsilon.S<-qunif(lhs.sample$epsilon.S,1.5,10)
lhs.sample$Gamma_re<-qunif(lhs.sample$Gamma_re,0.001,0.2)
lhs.sample$Gamma_ce<-qunif(lhs.sample$Gamma_ce,0.001,0.2)
lhs.sample$Eta_re.N<-qunif(lhs.sample$Eta_re.N,0.5,0.95)
lhs.sample$Eta_ce.N<-qunif(lhs.sample$Eta_ce.N,0.5,0.95)
lhs.sample$Eta_re.S<-qunif(lhs.sample$Eta_re.S,0.001,0.10)
lhs.sample$Eta_ce.S<-qunif(lhs.sample$Eta_ce.S,0.001,0.10)

# Para val.param.N y val.param.S, asumo entorno ~1.0 ±50%
lhs.sample$val.param.N <- qunif(lhs.sample$val.param.N, 0.5, 1.5)
lhs.sample$val.param.S <- qunif(lhs.sample$val.param.S, 0.5, 1.5)

# Para alfa.N y alfa.S, asumo entre 0.1 y 0.9
lhs.sample$alfa.N     <- qunif(lhs.sample$alfa.N,     0.1, 0.9)
lhs.sample$alfa.S     <- qunif(lhs.sample$alfa.S,     0.1, 0.9)
#

#multi-core runs
library(snow)
library(mFilter)       # para que tu sesión principal conozca hpfilter()
#  nCore<-40
 nCore<-8
 cl <- makeSOCKcluster(names = rep('localhost',nCore))
 
 # 1) Carga mFilter en cada nodo, para que hpfilter() exista allí:
 clusterEvalQ(cl, {
   library(mFilter)
   NULL
 })
 
 global.elements<-list("ediamMSPE","ediamCalib","ediamEmpiricalParameters","ediamInitialConditions","EdiamEquations","Ediam","HistData","ValVars","dir.model")
 clusterExport(cl,global.elements,envir=environment())

 test.out<- parApply(cl,lhs.sample,1,function(x) {out<-ediamMSPE(c( as.numeric(x['epsilon.N']),
                                                                   as.numeric(x['epsilon.S']),
                                                                   as.numeric(x['Gamma_re']),
                                                                   as.numeric(x['Gamma_ce']),
                                                                   as.numeric(x['Eta_re.N']),
                                                                   as.numeric(x['Eta_ce.N']),
                                                                   as.numeric(x['Eta_re.S']),
                                                                   as.numeric(x['Eta_ce.S']),
                                                                   as.numeric(x['val.param.N']),
                                                                   as.numeric(x['val.param.S']),
                                                                   as.numeric(x['alfa.N']),
                                                                   as.numeric(x['alfa.S'])
                                                                 ),verbose=TRUE);
                                                 out$Run.ID<-as.numeric(x['Run.ID'])
                                       return(out)})
#
 stopCluster(cl)

 test.out<-do.call("rbind",test.out)

#merge with parametrs
dim(test.out)
dim(lhs.sample)
 test.out<-merge(test.out,lhs.sample, by = "Run.ID")
dim(test.out)

test.out<-test.out[order(test.out$MSPE.all ),]
head(test.out)





  
  

  
  








