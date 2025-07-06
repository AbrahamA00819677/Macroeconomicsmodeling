#in pc. Este trabajo esta hecho sobre la base del trabajo del Dr. Edmundo Molina-Perez
root<-"~/Documents/macroeconomia/proyectoFinal"
model.version<-"Ediam_v2020_02_18"
dir.calib <- file.path(root, model.version, "CalibrationScripts")
dir.model <- file.path(root, model.version)
dir.data  <- file.path(dir.calib, "CalibrationData")


#=============================================================================================
# Experiment I: Using historical data from 1995 to 2014 (Calibration)
#============================================================================================
#Load library for using HP filter
library(mFilter)
#
#Load script for running the historical calibration
source( file.path(dir.calib, "ediam_HistoricCalib.r"))
#Load script for determining initial conditions
source( file.path(dir.calib, "ediam_InitialConditions.r") )
#Load EDIAM Equations
source( file.path(dir.model, "ediam_Equations_discrete.r") )

#load historical data
HistData <<- read.csv( file.path(dir.data, "HistData.csv") )
HistData<<-subset(HistData,HistData$Time%in%c(1971:2014))
#Subset the data to the years that we have data for
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
#HistData$c_oil<-HistData[,"OILTPES"]/(HistData[,"TPES"]*HistData[,"FOSSILTPES"])#[1]
HistData$c_oil<-HistData[,"OILTPES"]/(HistData[,"Y_ce.h"])#[1]
HistData$Roil.h<-HistData[,"OILTPES"] #MTOE
HistData$Re.h<-HistData$Roil.h/HistData$c_oil #MTOE #this is something else I could use to calibrate, we could compare model's output against this time series
HistData$Price.Oil<- HistData$Price.Oil*7.1428571428571*1e6/1e9 # (billion USD per Mtoe) ; Asumming that 1 toe= 7.1428571428571 boe
HistData$ReToGDP.h<-(HistData$Re.h*HistData$Price.Oil)/(HistData$GDP) # [1] #are these two (Price Oil and GDP) in real terms ????
#hp filter to oil prices
HistData$Price.Oil.hp<-rep(as.numeric(hpfilter(subset(HistData[,"Price.Oil"],HistData[,"Country"]=="World"),freq=100)$trend),length(unique(HistData$Country)))


#Indicate variables of interest for estimating MSPE
ValVars<<-c("GDP.N","Y_ce.N","Y_re.N","GDP.S","Y_ce.S","Y_re.S")

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
    lrng.re.N =  round(x[11],3), # #Note: no experience accumulation effects considered in this experiment
    lrng.ce.N = round(x[12],3), #0.0, #Note: no experience accumulation effects considered in this experiment
    lrng.re.S = round(x[13],3), #0.0, #Note: no experience accumulation effects considered in this experiment
    lrng.ce.S = round(x[14],3), #0.0, #Note: no experience accumulation effects considered in this experiment
    alfa.N = round(x[15],4),
    alfa.S = round(x[16],4)
  )
  
  #historic data with hp
  Y_re.Nh <<-as.numeric(hpfilter(subset(HistData[,"Y_re.h"],HistData[,"Country"]=="Memo: OECD Total"),freq=100)$trend) #Mtoe
  Y_ce.Nh <<-as.numeric(hpfilter(subset(HistData[,"Y_ce.h"],HistData[,"Country"]=="Memo: OECD Total"),freq=100)$trend) #Mtoe
  Y_re.Sh <<-as.numeric(hpfilter(subset(HistData[,"Y_re.h"],HistData[,"Country"]=="NonOECD"),freq=100)$trend) #Mtoe
  Y_ce.Sh <<-as.numeric(hpfilter(subset(HistData[,"Y_ce.h"],HistData[,"Country"]=="NonOECD"),freq=100)$trend) #Mtoe
  
  #Oil prices
  Price.oil.y<<- subset(HistData[,"Price.Oil.hp"],HistData[,"Country"]=="World")
  
  #GDP with hp
  GDP.Nh<<-as.numeric(hpfilter(subset(HistData[,"GDP"],HistData[,"Country"]=="Memo: OECD Total"),freq=100)$trend) #(billion 2010 USD using exchange rates)
  GDP.Sh<<-as.numeric(hpfilter(subset(HistData[,"GDP"],HistData[,"Country"]=="NonOECD"),freq=100)$trend) #(billion 2010 USD using exchange rates)
  
  #Oil Supply
  Re.Nh<<-as.numeric(hpfilter(subset(HistData[,"Re.h"],HistData[,"Country"]=="Memo: OECD Total"),freq=100)$trend) # Mtoe
  Re.Sh<<-as.numeric(hpfilter(subset(HistData[,"Re.h"],HistData[,"Country"]=="NonOECD"),freq=100)$trend) #Mtoe
  ReToGDP.Nh<<-(Re.Nh*Price.oil.y)/(GDP.Nh) # [1]
  ReToGDP.Sh<<-(Re.Sh*Price.oil.y)/(GDP.Sh) # [1]
  
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
  #   tolerance<-0.11
  #   MSPE<- ifelse(mean(MSPEp1)<=tolerance,
  #               ifelse(mean(MSPEp1.1)<=tolerance,
  #                       ifelse(mean(MSPEp1.2)<=tolerance,
  #                                ifelse(mean(MSPEp2)<=tolerance,
  #                                     ifelse(mean(MSPEp2.1)<=tolerance,
  #                                          ifelse(mean(MSPEp2.2)<=tolerance,MSPE,1e3)
  #                                                                 ,1e3)
  #                                                        ,1e3)
  #                                                    ,1e3)
  #                                            ,1e3)
  #                                    ,1e3)
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

# Use algorithm genoud to estimate parameters
#############################################################################################
#Optimization
#Load library for parallel computing
library(snow)
#Specify numbers of cores available for calibration
nCore<- parallel::detectCores() - 1 #about 50 cores is the optimal and it takes an hour to run, but we could run severl splits of 50,
#Define cluster
cl <- makeSOCKcluster(names = rep('localhost',nCore))
global.elements<-list("ediamMSPE","ediamCalib","ediamEmpiricalParameters","ediamInitialConditions","EdiamEquations","Ediam","HistData","ValVars","dir.model","hpfilter")
clusterExport(cl,global.elements,envir=environment())

#Load libary genoud
library(rgenoud)
#set seed for genetic optimization
set.seed(55555)
#Execute the optimization
out<-genoud(ediamMSPE,max=FALSE,
            nvars=16,
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
                                 0.000 , #lrng.re.N
                                 0.000 ,#, #lrng.ce.N
                                 0.000 , #lrng.re.N
                                 0.000 , #lrng.ce.S
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
              0.000 , #lrng.re.N
              0.000 , #lrng.ce.N
              0.000 , #lrng.re.S
              0.000 , #lrng.ce.S,
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
              1.0 , #lrng.re.N
              1.0 ,#, #lrng.ce.N
              1.0 , #lrng.re.S
              1.0 , #lrng.ce.S
              0.5, #alfa.N
              0.5 #alfa.S
            ),
            ncol=2),
            boundary.enforcement = 1,
            cluster=cl,
            print.level=1)

stopCluster(cl)

print(out$par)
print(out$value)
#out

errors<- ediamMSPE(out$par, verbose = TRUE)
#save calibrated parameters
calib.params<<-c(
  epsilon.N = round(out$par[1],4),
  epsilon.S = round(out$par[2],4),
  Gamma_re =  round(out$par[3],4),
  Gamma_ce=   round(out$par[4],4),
  Eta_re.N =  round(out$par[5],4),
  Eta_ce.N =  round(out$par[6],4),
  Eta_re.S =  round(out$par[7],4),
  Eta_ce.S =  round(out$par[8],4),
  val.param.N = round(out$par[9],4),
  val.param.S = round(out$par[10],4),
  lrng.re.N = round(out$par[11],4), 
  lrng.ce.N = round(out$par[12],4), 
  lrng.re.S = round(out$par[13],4), 
  lrng.ce.S = round(out$par[14],4), 
  alfa.N = round(out$par[15],4),
  alfa.S = round(out$par[16],4)
)

#Create historical reference
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


dir.data <- file.path(root, model.version, "Data")

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
  
  # Vector en entorno global
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
global.elements<-list("ediamPolicyEval","ediamMain","ediamEmpiricalParameters","ediamInitialConditions","EdiamEquations","Ediam","oil.scenario","pop.scenario","climate.scenario","calib.params" )
clusterExport(cl,global.elements,envir=environment())

#Load libary genoud
library(rgenoud)
#set seed for genetic optimization
set.seed(55555)
#Execute the optimization
out2<-genoud(ediamPolicyEval,max=TRUE,
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
print(out2$par)
print(out2$value)

x  <- out2$par       # tu vector de 12 policy-parameters

SimulData<-ediamPolicyEval(x,verbose=TRUE)

library(ggplot2)
library(reshape2)
library(ggpubr)

ediamPlot(SimulData)


#This function plots results from run in Ediam

ediamPlot4<-function(SimulData)
{
  #create function to use data by region
  data.by.region<-function(data,id.vars,measure.vars)
  {
    data<-data[,c(id.vars,measure.vars)]
    data<-melt(data, id.vars=id.vars, measure.vars=measure.vars,  variable.name="Variable.Region")
    data<-cbind(data,colsplit(data$Variable.Region,("\\."),c("Variable","Region")))
    data$Variable.Region<-NULL
    new.form<-as.formula(paste(id.vars,"Region~Variable",sep="+"))
    data<-dcast(data,new.form, value.var="value")
    data$Region<-ifelse(data$Region=="N","OECD Countries","Non-OECD Countries")
    return(data)
  }
  
  
  #
  theme_paper<- theme( panel.border = element_rect(fill = NA, colour = "black", size = 1),
                       #text=element_text(family="Arial"),
                       #style of title
                       plot.title = element_text(size=11),
                       #style axis titles
                       axis.title.x = element_text(size=10),
                       axis.title.y = element_text(size=10),
                       #style axis
                       axis.text.x = element_text(size=8),
                       axis.text.y = element_text(size=8),
                       panel.grid.major = element_blank(),
                       panel.grid.major.x = element_blank(),
                       panel.grid.minor = element_blank(),
                       panel.grid.minor.x = element_blank(),
                       panel.background = element_blank())
  #
  Plot0 <- ggplot(SimulData, aes(times+2014,Delta.Temp))+ geom_line()+geom_line(size=1)+
    geom_hline(yintercept = 2.0,linetype = 2)+geom_hline(yintercept = 6.0,linetype = 2)+
    ggtitle("Global Temperature Anomaly With Respect to Historic Mean")+
    scale_x_continuous(name="Time [Years]", breaks=seq(2014,2014+100,by=10)) +
    scale_y_continuous(name="Degrees Celsius",limits=c(0, 6.0))+theme_paper
  #
  Plot1<-ggplot(data.by.region(SimulData,"times",c("Realcetax.N","Realcetax.S")), aes(times+2014,Realcetax,colour =Region)) + geom_line()+geom_line(size=1)+
    ggtitle("Carbon Price")+
    #scale_color_manual(values = c("deepskyblue3","goldenrod2"))+
    scale_color_manual(values =c("#E7B800","#00AFBB"))+
    scale_x_continuous(name="Time [Years]", breaks=seq(2014,2014+100,by=10)) +
    scale_y_continuous(name="USD 2010 per tonne of co2")+theme_paper
  #
  Plot2<-ggplot(data.by.region(SimulData,"times",c("RealRDsubsidy.N","RealRDsubsidy.S")), aes(times+2014,RealRDsubsidy,colour =Region)) + geom_line()+geom_line(size=1)+
    ggtitle("Renewable Energy R&D Investments")+
    #scale_color_manual(values = c("deepskyblue3","goldenrod2"))+
    scale_color_manual(values =c("#E7B800","#00AFBB"))+
    scale_x_continuous(name="Time [Years]", breaks=seq(2014,2014+100,by=10)) +
    scale_y_continuous(name="Billion USD 2010")+theme_paper
  
  #
  Plot3<-ggplot(data.by.region(SimulData,"times",c("RealTecsubsidy.N","RealTecsubsidy.S")), aes(times+2014,RealTecsubsidy,colour =Region)) + geom_line()+geom_line(size=1)+
    ggtitle("Renewable Energy Investments")+
    #scale_color_manual(values = c("deepskyblue3","goldenrod2"))+
    scale_color_manual(values =c("#E7B800","#00AFBB"))+
    scale_x_continuous(name="Time [Years]", breaks=seq(2014,2014+100,by=10)) +
    scale_y_continuous(name="Billion USD 2010")+theme_paper
  
  #
  #Graph Scientists
  Plot4<-ggplot(data.by.region(SimulData,"times",c("s_re.N","s_re.S")), aes(times+2014,s_re,colour =Region)) + geom_line()+geom_line(size=1)+
    ggtitle("Entrepreneurs in Renewable Energy")+
    #scale_color_manual(values = c("deepskyblue3","goldenrod2"))+
    scale_color_manual(values =c("#E7B800","#00AFBB"))+
    scale_x_continuous(name="Time [Years]", breaks=seq(2014,2014+100,by=10)) +
    scale_y_continuous(name="Percent of Total", breaks=seq(0,1,by=0.20),limits=c(0, 1.0))+theme_paper
  
  #
  #Graph Rate of decarbonization
  Plot5<-ggplot(data.by.region(SimulData,"times",c("DecarbY.N","DecarbY.S")), aes(times+2014,DecarbY,colour =Region)) + geom_line(size=1)+
    ggtitle("Decarbonization Rate")+
    #scale_color_manual(values = c("deepskyblue3","goldenrod2"))+
    scale_color_manual(values =c("#E7B800","#00AFBB"))+
    scale_x_continuous(name="Time [Years]", breaks=seq(2014,2014+100,by=10)) +
    scale_y_continuous(name="Percent of Total", breaks=seq(0,1,by=0.20),limits=c(0, 1.0))+theme_paper
  
  
  #
  #World Population
  Plot6<-ggplot(data.by.region(SimulData,"times",c("L.N","L.S")), aes(times+2014,L,fill =Region)) +geom_bar(stat="identity")+
    ggtitle("Population Growth")+
    scale_fill_manual(values =c("#E7B800","#00AFBB"))+
    scale_x_continuous(name="Time [Years]", breaks=seq(2014,2014+100,by=10)) +
    scale_y_continuous(name="Million People")+theme_paper
  
  #Price of Renewable Energy
  Plot7<-ggplot(data.by.region(SimulData,"times",c("Price_re.N","Price_re.S")), aes(times+2014,Price_re,fill =Region)) +geom_bar(stat="identity")+
    ggtitle("Price of Renewable Energy")+
    scale_fill_manual(values =c("#E7B800","#00AFBB"))+
    scale_x_continuous(name="Time [Years]", breaks=seq(2014,2014+100,by=10)) +
    scale_y_continuous(name="USD")+theme_paper
  #Price of Renewable Energy
  Plot8<-ggplot(data.by.region(SimulData,"times",c("Profits_re.N","Profits_re.S")), aes(times+2014,Profits_re,fill =Region)) +geom_bar(stat="identity")+
    ggtitle("Profits of Renewables productors")+
    scale_fill_manual(values =c("#E7B800","#00AFBB"))+
    scale_x_continuous(name="Time [Years]", breaks=seq(2014,2014+100,by=10)) +
    scale_y_continuous(name="USD")+theme_paper
  #Production of Renewable Energy
  Plot9<-ggplot(data.by.region(SimulData,"times",c("Y_re.N","Y_re.S")), aes(times+2014,Y_re,fill =Region)) +geom_bar(stat="identity")+
    ggtitle("Production of Renewable Energy")+
    scale_fill_manual(values =c("#E7B800","#00AFBB"))+
    scale_x_continuous(name="Time [Years]", breaks=seq(2014,2014+100,by=10)) +
    scale_y_continuous(name="energy")+theme_paper
  #Price of fossil energy
  Plot10<-ggplot(data.by.region(SimulData,"times",c("Price_ce.N","Price_ce.S")), aes(times+2014,Price_ce,fill =Region)) +geom_bar(stat="identity")+
    ggtitle("Price of fossil energy")+
    scale_fill_manual(values =c("#E7B800","#00AFBB"))+
    scale_x_continuous(name="Time [Years]", breaks=seq(2014,2014+100,by=10)) +
    scale_y_continuous(name="USD")+theme_paper
  #Production of Renewable Energy
  Plot11<-ggplot(data.by.region(SimulData,"times",c("Profits_ce.N","Profits_ce.S")), aes(times+2014,Profits_ce,fill =Region)) +geom_bar(stat="identity")+
    ggtitle("Profits of fossil productors")+
    scale_fill_manual(values =c("#E7B800","#00AFBB"))+
    scale_x_continuous(name="Time [Years]", breaks=seq(2014,2014+100,by=10)) +
    scale_y_continuous(name="USD")+theme_paper
  
  #Production of fossil Energy
  Plot12<-ggplot(data.by.region(SimulData,"times",c("Y_ce.N","Y_ce.S")), aes(times+2014,Y_ce,fill =Region)) +geom_bar(stat="identity")+
    ggtitle("Production of fossil Energy")+
    scale_fill_manual(values =c("#E7B800","#00AFBB"))+
    scale_x_continuous(name="Time [Years]", breaks=seq(2014,2014+100,by=10)) +
    scale_y_continuous(name="energy")+theme_paper
  
  #Production of total Energy
  Plot13<-ggplot(data.by.region(SimulData,"times",c("Y.N","Y.S")), aes(times+2014,Y,fill =Region)) +geom_bar(stat="identity")+
    ggtitle("Production of total Energy")+
    scale_fill_manual(values =c("#E7B800","#00AFBB"))+
    scale_x_continuous(name="Time [Years]", breaks=seq(2014,2014+100,by=10)) +
    scale_y_continuous(name="energy")+theme_paper
  
  
  #Price of total energy
  Plot14<-ggplot(data.by.region(SimulData,"times",c("Price_Y.N","Price_Y.S")), aes(times+2014,Price_Y,fill =Region)) +geom_bar(stat="identity")+
    ggtitle("Price of total energy")+
    scale_fill_manual(values =c("#E7B800","#00AFBB"))+
    scale_x_continuous(name="Time [Years]", breaks=seq(2014,2014+100,by=10)) +
    scale_y_continuous(name="USD")+theme_paper
  
  #GDP
  Plot15<-ggplot(data.by.region(SimulData,"times",c("GDP.N","GDP.S")), aes(times+2014,GDP,fill =Region)) +geom_bar(stat="identity")+
    ggtitle("GDP")+
    scale_fill_manual(values =c("#E7B800","#00AFBB"))+
    scale_x_continuous(name="Time [Years]", breaks=seq(2014,2014+100,by=10)) +
    scale_y_continuous(name="USD")+theme_paper
  
  #GrowthRate
  Plot16<-ggplot(data.by.region(SimulData,"times",c("GrowthRate.N","GrowthRate.S")), aes(times+2014,GrowthRate,fill =Region)) +geom_bar(stat="identity")+
    ggtitle("GrowthRate.S")+
    scale_fill_manual(values =c("#E7B800","#00AFBB"))+
    scale_x_continuous(name="Time [Years]", breaks=seq(2014,2014+100,by=10)) +
    scale_y_continuous(name="USD")+theme_paper
  
  #Consumption
  Plot17<-ggplot(data.by.region(SimulData,"times",c("Consumption.N","Consumption.S")), aes(times+2014,Consumption,fill =Region)) +geom_bar(stat="identity")+
    ggtitle("Consumption aggregate")+
    scale_fill_manual(values =c("#E7B800","#00AFBB"))+
    scale_x_continuous(name="Time [Years]", breaks=seq(2014,2014+100,by=10)) +
    scale_y_continuous(name="USD")+theme_paper
  
  
  
  
  
  #Putting all graphs into one page
  
  Plot<-ggarrange(Plot1,Plot4,
                  Plot2,Plot5,
                  Plot3,Plot6,Plot7,Plot8,Plot9,Plot10,Plot11,Plot12,Plot13,Plot14, Plot15,Plot16, Plot17,
                  Plot0,
                  nrow=9,ncol = 2,
                  common.legend = TRUE,
                  align="h")
  
  return(Plot)
}

ediamPlot4(SimulData)



## Experimentacion, transicion energetica enfocada en la minimiuzacion del delta.temp

## Definimos nueva funcion para optimizar el delta de temperatura

ediamPolicyTempEval <- function(x, verbose = FALSE) {
  
  # Paso de simulación
  TimeStep <<- 5
  
  # Datos históricos y escenarios
  Y_re.Nh <<- 490.0097;  Y_ce.Nh <<- 4255.120
  Y_re.Sh <<- 1386.1706; Y_ce.Sh <<- 7002.889
  Re.Nh    <<- 4255.120;  Re.Sh    <<- 7002.889
  ReToGDP.Nh <<- 0.06227680; ReToGDP.Sh <<- 0.18765958
  
  oil.times   <<- oil.scenario$Year - 2014
  Price.oil.y <<- oil.scenario$Price.Oil
  
  GDP.Nh <<- 47034.07
  GDP.Sh <<- 25688.192
  
  pop.times <<- pop.scenario$Year - 2014
  L.N.y     <<- pop.scenario$AdvancedRegion
  L.S.y     <<- pop.scenario$EmergingRegion
  
  # Parámetros climáticos
  climate.params <<- c(
    qsi                  = climate.scenario$qsi,
    Delta.S              = climate.scenario$Delta.S,
    Delta.Temp.Disaster  = climate.scenario$Delta.Temp.Disaster,
    Beta.Delta.Temp      = climate.scenario$Beta.Delta.Temp,
    CO2.base             = climate.scenario$CO2.base,
    CO2.Disaster         = climate.scenario$CO2.Disaster,
    DeltaTempConstant    = climate.scenario$DeltaTempConstant,
    S_0                  = climate.scenario$S.0
  )
  
  # Vector de políticas (12 parámetros)
  policy.params <<- c(
    ce.tax.N           = round(x[1],  3),
    Schedule.ce.tax.N  = round(x[2],  4),
    Tec.subsidy.N      = round(x[3],  3),
    Schedule.Tec.subsidy.N = round(x[4],  4),
    RD.subsidy.N       = round(x[5],  3),
    Schedule.RD.subsidy.N  = round(x[6],  4),
    ce.tax.S           = round(x[7],  3),
    Schedule.ce.tax.S  = round(x[8],  4),
    Tec.subsidy.S      = round(x[9],  3),
    Schedule.Tec.subsidy.S = round(x[10], 4),
    RD.subsidy.S       = round(x[11], 3),
    Schedule.RD.subsidy.S  = round(x[12], 4)
  )
  
  #–– 2) Ejecución y captura de fallos ––
  result <- tryCatch({
    sim <- ediamMain(calib.params, verbose = TRUE)
    finalTemp <- tail(sim$Delta.Temp, 1)
    
    if (!is.finite(finalTemp) || is.na(finalTemp)) {
      # Penalizamos valores inválidos
      return(1e6)
    }
    
    if (!verbose) {
      return(finalTemp)
    } else {
      return(sim)
    }
    
  }, error = function(e) {
    # Si hay cualquier error, devolvemos un gran costo
    return(1e6)
  })
  
  return(result)
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
global.elements<-list("ediamPolicyEval","ediamMain","ediamEmpiricalParameters","ediamInitialConditions","EdiamEquations","Ediam","oil.scenario","pop.scenario","climate.scenario","calib.params" )
clusterExport(cl,global.elements,envir=environment())

#Load libary genoud
library(rgenoud)
#set seed for genetic optimization
set.seed(55555)
#Execute the optimization
out3<-genoud(ediamPolicyTempEval,max=FALSE,
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
             optim.method  = "Nelder-Mead",  # desactiva BFGS
             gradient.check = FALSE,
             print.level=1)

stopCluster(cl)
print(out3$par)
print(out3$value)

x2  <- out3$par       # tu vector de 12 policy-parameters

SimulData<-ediamPolicyEval(x2,verbose=TRUE)


ediamPlot4(SimulData)

### Tercer experimento equiliobrio entre welfare y DeltaTemp

# Instala e importa mco si no lo tienes
# install.packages("mco")
library(mco)

# Función que devuelve un vector de dos objetivos: 
# (1) -Welfare  y (2) Temperatura final
ediamMultiObj <- function(x) {
  # Configura exactamente igual que en ediamPolicyTempEval / original
  TimeStep <<- 5
  oil.times   <<- oil.scenario$Year - 2014
  Price.oil.y <<- oil.scenario$Price.Oil
  pop.times   <<- pop.scenario$Year - 2014
  L.N.y       <<- pop.scenario$AdvancedRegion
  L.S.y       <<- pop.scenario$EmergingRegion
  climate.params <<- c(
    qsi                 = climate.scenario$qsi,
    Delta.S             = climate.scenario$Delta.S,
    Delta.Temp.Disaster = climate.scenario$Delta.Temp.Disaster,
    Beta.Delta.Temp     = climate.scenario$Beta.Delta.Temp,
    CO2.base            = climate.scenario$CO2.base,
    CO2.Disaster        = climate.scenario$CO2.Disaster,
    DeltaTempConstant   = climate.scenario$DeltaTempConstant,
    S_0                 = climate.scenario$S.0
  )
  policy.params <<- c(
    ce.tax.N            = round(x[1],  3),
    Schedule.ce.tax.N   = round(x[2],  4),
    Tec.subsidy.N       = round(x[3],  3),
    Schedule.Tec.subsidy.N = round(x[4], 4),
    RD.subsidy.N        = round(x[5],  3),
    Schedule.RD.subsidy.N  = round(x[6], 4),
    ce.tax.S            = round(x[7],  3),
    Schedule.ce.tax.S   = round(x[8],  4),
    Tec.subsidy.S       = round(x[9],  3),
    Schedule.Tec.subsidy.S = round(x[10],4),
    RD.subsidy.S        = round(x[11], 3),
    Schedule.RD.subsidy.S  = round(x[12],4)
  )
  
  # Corre la simulación en modo verbose para extraer todo
  sim <- ediamMain(calib.params, verbose = TRUE)
  
  # Objetivo 1: NEGATIVO de welfare (pues nsga2 minimiza)
  W    <- sim$Social.Welfare.Function[1]
  f1   <- -W
  
  # Objetivo 2: Temperatura final
  Tfin <- tail(sim$Delta.Temp, 1)
  f2   <- Tfin
  
  return(c(f1, f2))
}

# Parámetros para nsga2
nvars       <- 12
lower.bounds <- rep(0, nvars)
upper.bounds <- c(1.0,0.1,0.9,0.1,4.0,0.1, 1.0,0.1,0.9,0.1,4.0,0.1)

# Ejecuta NSGA‐II
set.seed(12345)
res <- nsga2(
  fn        = ediamMultiObj, 
  idim      = nvars, 
  odim      = 2,             # dos objetivos
  lower.bounds = lower.bounds, 
  upper.bounds = upper.bounds,
  popsize   = 60,           # tamaño de población
  generations = 20          # número de generaciones
)

# Extraer Pareto front
pareto.objs <- res$value       # matriz N×2 con (f1, f2)
pareto.pars <- res$par         # matriz N×12 con los x asociados

# Reconviértelos a (W, T) para interpretación
welfare   <- -pareto.objs[,1]
temp.final <- pareto.objs[,2]

# Visualiza el frente
plot(welfare, temp.final, 
     xlab = "Welfare", ylab = "Δ Temp final",
     main = "Fronte de Pareto: Welfare vs Temperatura")

# Opcional: selecciona un punto del frente (por ejemplo el que maximiza Welfare - Temperatura)
scores <- welfare - temp.final
best.idx <- which.max(scores)
best.x   <- pareto.pars[best.idx, ]



#Misma funcion que la de arriba, ase agrega el procesamiento en paralelo para eficientizar la busqueda
# ================================
# 1) Carga de librerías
# ================================
if (!requireNamespace("mco", quietly = TRUE)) install.packages("mco")
library(parallel)
library(snow)
library(mco)

# ================================
# 2) Función multi‐objetivo autocontenida
#    devuelve c(-W, ΔTemp_final)
# ================================
ediamMultiObj <- function(x) {
  # —– 2a) Datos históricos y step ––
  TimeStep <<- 5
  Y_re.Nh    <<- 490.0097;  Y_ce.Nh    <<- 4255.120
  Y_re.Sh    <<- 1386.1706; Y_ce.Sh    <<- 7002.889
  Re.Nh      <<- 4255.120;  Re.Sh      <<- 7002.889
  ReToGDP.Nh <<- 0.06227680; ReToGDP.Sh <<- 0.18765958
  GDP.Nh     <<- 47034.07
  GDP.Sh     <<- 25688.192
  
  # —– 2b) Escenarios exógenos ––
  oil.times   <<- oil.scenario$Year - 2014
  Price.oil.y <<- oil.scenario$Price.Oil
  pop.times   <<- pop.scenario$Year - 2014
  L.N.y       <<- pop.scenario$AdvancedRegion
  L.S.y       <<- pop.scenario$EmergingRegion
  
  # —– 2c) Parámetros climáticos ––
  climate.params <<- c(
    qsi                 = climate.scenario$qsi,
    Delta.S             = climate.scenario$Delta.S,
    Delta.Temp.Disaster = climate.scenario$Delta.Temp.Disaster,
    Beta.Delta.Temp     = climate.scenario$Beta.Delta.Temp,
    CO2.base            = climate.scenario$CO2.base,
    CO2.Disaster        = climate.scenario$CO2.Disaster,
    DeltaTempConstant   = climate.scenario$DeltaTempConstant,
    S_0                 = climate.scenario$S.0
  )
  
  # —– 2d) Vector de políticas ––
  policy.params <<- c(
    ce.tax.N              = round(x[1],  3),
    Schedule.ce.tax.N     = round(x[2],  4),
    Tec.subsidy.N         = round(x[3],  3),
    Schedule.Tec.subsidy.N= round(x[4],  4),
    RD.subsidy.N          = round(x[5],  3),
    Schedule.RD.subsidy.N = round(x[6],  4),
    ce.tax.S              = round(x[7],  3),
    Schedule.ce.tax.S     = round(x[8],  4),
    Tec.subsidy.S         = round(x[9],  3),
    Schedule.Tec.subsidy.S= round(x[10], 4),
    RD.subsidy.S          = round(x[11], 3),
    Schedule.RD.subsidy.S = round(x[12], 4)
  )
  
  # —– 2e) Corre la simulación ––
  sim <- ediamMain(calib.params, verbose = TRUE)
  
  # —– 2f) Objetivos ––
  W    <- sim$Social.Welfare.Function[1] # welfare
  f1   <- -W                             # lo minimiza nsga2 => maximiza W
  f2   <- tail(sim$Delta.Temp, 1)       # ΔTemp final
  
  return(c(f1, f2))
}

# ================================
# 3) Montaje del clúster snow
# ================================
nCore <- detectCores() - 1
cl    <- makeSOCKcluster(rep("localhost", nCore))

clusterExport(cl,
              varlist = c(
                "ediamMultiObj",
                "ediamMain","ediamEmpiricalParameters","ediamInitialConditions","EdiamEquations","Ediam","oil.scenario","pop.scenario","climate.scenario","calib.params"
              ),
              envir = environment()
)
clusterEvalQ(cl, library(mco))  # si ediamMain usa mco internamente

# ================================
# 4) Wrapper vectorizado para nsga2
# ================================
vectorizedEval <- function(pop) {
  # parApply(…) ya es 2×popsize
  parApply(cl, pop, 1, ediamMultiObj)
}


# ================================
# 5) Ejecuta NSGA‐II
# ================================
set.seed(12345)
res <- nsga2(
  fn           = vectorizedEval,
  idim         = 12,
  odim         = 2,
  lower.bounds = rep(0, 12),
  upper.bounds = c(1.0,0.1,0.9,0.1,4.0,0.1,
                   1.0,0.1,0.9,0.1,4.0,0.1),
  popsize      = 100,
  generations  = 50,
  vectorized   = TRUE
)

# ================================
# 6) Limpieza y visualización
# ================================
stopCluster(cl)

pareto.objs <- res$value
pareto.pars <- res$par

welfare    <- -pareto.objs[,1]
temp.final <-  pareto.objs[,2]

plot(welfare, temp.final,
     xlab="Welfare", ylab="Δ Temp final",
     main="Fronte de Pareto: Welfare vs Temperatura", pch=19)




# 1) Encuentra los índices con valores válidos en ambos vectores
valid_idx <- which(is.finite(welfare) & is.finite(temp.final))

# 2) Si no hay ninguno, aborta
if (length(valid_idx) == 0) {
  stop("No hay puntos válidos en el frente de Pareto.")
}

# 3) Punto ideal sobre ese subconjunto
idealW <- max(welfare[valid_idx])
idealT <- min(temp.final[valid_idx])

# 4) Distancias sólo para los índices válidos
dists2 <- sqrt(
  (welfare[valid_idx] - idealW)^2 +
    (temp.final[valid_idx] - idealT)^2
)

# 5) Índice “local” (en valid_idx) del mejor
which_loc   <- which.min(dists2)

# 6) Convertirlo al índice “global” sobre todo el frente
best_global <- valid_idx[which_loc]

# 7) Imprime resultados
cat(
  "Índice del mejor compromiso (global):", best_global, "\n",
  "  Welfare =", welfare[best_global], "\n",
  "  ΔT      =", temp.final[best_global], "\n"
)

# 8) También puedes ver la política asociada:
best_policy <- pareto.pars[best_global, ]
print(best_policy)



SimulData_opt <- ediamPolicyTempEval(best_policy, verbose = TRUE)
ediamPlot4(SimulData_opt)


