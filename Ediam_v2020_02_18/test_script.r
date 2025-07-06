#we need to run a test, use the model to run a simulation with two regions with the same parameters, make sure there isn't an error in the model
#review how the stock of labour is computed, dicrepancies in accounting for the size of the economy might be creating a problem

#in pc
  root<-"C:\\Users\\L03054557\\OneDrive\\Edmundo-ITESM\\3.Proyectos\\Disertation\\Model Disertation\\TechChange Model\\"
  model.version<-"Ediam_v2019_12_20"
  dir.calib <<- paste(root,model.version,"\\CalibrationScripts\\",sep="")
  dir.model <<- paste(root,model.version,"\\",sep="")
  dir.data <<- paste(root,model.version,"\\CalibrationScripts\\CalibrationData\\",sep="")


#load Model
source(paste(dir.calib,"ediam_HistoricCalib.r",sep=""))
#Load script for determining initial conditions
source(paste(dir.calib,"ediam_InitialConditions.r",sep=""))
#Load EDIAM Equations
source(paste(dir.model,"ediam_Equations_discrete.r",sep=""))

 #load calibration params
   dir.bootstrap<-paste0(dir.calib,"\\Bootstrap_Out\\")
   filenames <- list.files(dir.bootstrap, pattern="*.csv", full.names=FALSE)
   params<-lapply(filenames, function (x) {read.csv(paste(dir.bootstrap,x,sep=""))})
   params<-do.call("rbind",params)



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

#two exact regions
 #Y_re.Nh <- Y_re.Sh
 #Y_ce.Nh <- Y_ce.Sh
 #Re.Nh <- Re.Sh
 #ReToGDP.Nh <- ReToGDP.Sh
 #GDP.Nh <- GDP.Sh
 #L.N.y <- L.S.y

  i<-1
  pars<-c(
             epsilon.N = round(as.numeric(params[1,"epsilon.N"]),4),
             epsilon.S = round(as.numeric(params[1,"epsilon.S"]),4),
             Gamma_re =  round(as.numeric(params[1,"Gamma_re"]),4),
             Gamma_ce =   round(as.numeric(params[1,"Gamma_ce"]),4),
             Eta_re.N =  round(as.numeric(params[1,"Eta_re.N"]),4),
             Eta_ce.N =  round(as.numeric(params[1,"Eta_ce.N"]),4),
             Eta_re.S =  round(as.numeric(params[1,"Eta_re.S"]),4),
             Eta_ce.S =  round(as.numeric(params[1,"Eta_ce.S"]),4),
             val.param.N = round(as.numeric(params[1,"val.param.N"]),4),
             val.param.S = round(as.numeric(params[1,"val.param.S"]),4),
             lrng.re.N = 0.0, #Note: no experience accumulation effects considered in this experiment
             lrng.ce.N = 0.0, #Note: no experience accumulation effects considered in this experiment
             lrng.re.S = 0.0, #Note: no experience accumulation effects considered in this experiment
             lrng.ce.S = 0.0, #Note: no experience accumulation effects considered in this experiment
             alfa.N = round(as.numeric(params[1,"alfa.N"]),4),
             alfa.S = round(as.numeric(params[1,"alfa.S"]),4)
           )
#
## =======================================================================================================================================================================================================
## Exploratory Dynamic Integreated Assessment Model
## =======================================================================================================================================================================================================
#ediamCalib<-function (calib.params, verbose=FALSE)
#{

calib.params<-pars

#==========================================================================================================================================================================================================
#Define policy parameters
#==========================================================================================================================================================================================================
# policy.vector<-c(
    #Carbon tax
#      ce.tax.N=0.0,
#      ce.tax.S=0.0,
    #Technology Policy in Advanced Region
#      Tec.subsidy.N = 0.0,
#      RD.subsidy.N = 0.0,
#      policy.half.life.N = 0.0,
    ##Technology Policy in Emerging Region
#      Tec.subsidy.S = 0.0,
#      RD.subsidy.S = 0.0,
#      policy.half.life.S = 0.0 )
#
policy.vector<-c(
   #Advanced Region
     ce.tax.N = 0.0,
      Schedule.ce.tax.N = 0.0,
     Tec.subsidy.N = 0.0,
      Schedule.Tec.subsidy.N = 0.0,
     RD.subsidy.N = 0.0,
      Schedule.RD.subsidy.N = 0.0,
      #policy.half.life.N = as.numeric(policy.params['policy.half.life.N']),
   #Emerging Region
     ce.tax.S= 0.0,
      Schedule.ce.tax.S = 0.0,
     Tec.subsidy.S = 0.0,
      Schedule.Tec.subsidy.S = 0.0,
     RD.subsidy.S = 0.0,
      Schedule.RD.subsidy.S = 0.0
     #policy.half.life.S = as.numeric(policy.params['policy.half.life.S'])
                )

#==========================================================================================================================================================================================================

#==========================================================================================================================================================================================================
#Define model parameters
#==========================================================================================================================================================================================================
 Parameters <- c(Run.ID = 1,
                 EndTime = (length(Price.oil.y)-1),
                 TimeStep = 1,
                #economic parameters advanced region
                 alfa.N = as.numeric(calib.params['alfa.N']),
                 epsilon.N = as.numeric(calib.params['epsilon.N']),
                 val.param.N = as.numeric(calib.params['val.param.N']),
                 lrng.re.N = as.numeric(calib.params['lrng.re.N']),
                 lrng.ce.N = as.numeric(calib.params['lrng.ce.N']),

                #economic parameters emerging region
                 alfa.S = as.numeric(calib.params['alfa.S']),
                 epsilon.S = as.numeric(calib.params['epsilon.S']),
                 val.param.S = as.numeric(calib.params['val.param.S']),
                 lrng.re.S = as.numeric(calib.params['lrng.re.S']),
                 lrng.ce.S = as.numeric(calib.params['lrng.ce.S']),

                #Technological paramters
                 Gamma_re = as.numeric(calib.params['Gamma_re']),
                 Gamma_ce = as.numeric(calib.params['Gamma_ce']),
                 Eta_re.N= as.numeric(calib.params['Eta_re.N']),
                 Eta_ce.N= as.numeric(calib.params['Eta_ce.N']),
                 Eta_re.S = as.numeric(calib.params['Eta_re.S']),
                 Eta_ce.S= as.numeric(calib.params['Eta_ce.S']),

                #Environmental parameters
                 earth.atm.weigth=5.3e18,#[kg]
                 molar.mass.CO2=44.0095, #[g/mol]
                 molar.mass.air=28.97, # [g/mol]
                 qsi= 0.000377221064255638,  #CO2 in ppm per 1 MTOE of FOSSILTPES
                 Delta.S = 0.00136321421736207,
				         Delta.Temp.Disaster = 6.0 ,
				         Beta.Delta.Temp = 4.12546004261621 ,
				         CO2.base = 276.80705087267,
                 CO2.Disaster = 1248.68034919411,
                 DeltaTempConstant =-0.215095719270412 ,
				         rho = 0.008,
				         lambda.S = 0.1443,
				         sigma.utility = 2.0,
                 S_0 = 1248.68034919411-315.2187 #S_0=CO2.Disaster-CO2[1971]
                )
#Add additional parameters
  Parameters<-c(Parameters,
                phi.N = (1-as.numeric(Parameters['alfa.N']))*(1-as.numeric(Parameters['epsilon.N'])),
                phi.S = (1-as.numeric(Parameters['alfa.S']))*(1-as.numeric(Parameters['epsilon.S'])),
                policy.vector)

#==========================================================================================================================================================================================================

#==========================================================================================================================================================================================================
#Define Initial Conditions
#==========================================================================================================================================================================================================

# Estimate empirical parameters
 EmpiricalParams<-apply(
                         data.frame(index= seq(1,1,by=1)),
                         1,function(x){ediamEmpiricalParameters(Parameters,as.numeric(x))})
 EmpiricalParams<-do.call("rbind",EmpiricalParams)
 EmpiricalParams<-apply(EmpiricalParams,2,mean)
# Add to parameters
 Parameters<-c(Parameters,EmpiricalParams)

#Use estimated parameters for setting initial conditions
 InitialConditions <-ediamInitialConditions(Parameters)

#Add additional parameters
 Parameters<-c(Parameters,InitialConditions$Parameters)

#Initial conditions for state variables
 InitialConditions <- InitialConditions$StateVariables
 GDP.Nh[1]
 GDP.Sh[1]

#==========================================================================================================================================================================================================
#Load EDIAM Difference Equations Structure and Execute Simulation
#==========================================================================================================================================================================================================

#Load Data Series
#for Oil
 Price.oil.Data <<- approxfun( x = seq(0,Parameters['EndTime'],1), #times
                       y = Price.oil.y, # Actual data points
                       method = "linear",
                       rule = 2)

#for Population
 L.N.Data<<-approxfun( x = seq(0,Parameters['EndTime'],1), #times
                       y = L.N.y, # Actual data points
                       method = "linear",
                       rule = 2)
#
 L.S.Data<<-approxfun( x = seq(0,Parameters['EndTime'],1), #times
                      y = L.S.y, # Actual data points
                      method = "linear",
                      rule = 2)
#Run the model
 out <- Ediam(end.time = Parameters['EndTime'] ,
              time.step = 1,
              y = InitialConditions,
              func = EdiamEquations,
              Parameters = Parameters)

#compare results
paste0("North ","model:", out$GDP.N[1], " vs ", "data:", GDP.Nh[1])
paste0("South ","model:", out$GDP.S[1], " vs ", "data:", GDP.Sh[1])

#==========================================================================================================================================================================================================
# Post-processing data and validation of constraints
#==========================================================================================================================================================================================================

#Define objetive function
# Social.Welfare.Function<-sum(as.numeric(out$Utility.Consumer_N))+sum(as.numeric(out$Utility.Consumer_S))

#Check constraints are valid
#Consumption constraint -In no year consumption can be negative-
# Social.Welfare.Function<-ifelse(min(out$Consumption.N)<0,-1000000,
#                                 ifelse(min(out$Consumption.S)<0,-1000000,Social.Welfare.Function))

#Budget constraint -the cost of policies cannot exceed total budget-
# Social.Welfare.Function<-ifelse(sum(out$Budget.function.N)<0,-1000000,
#                                 ifelse(sum(out$Budget.function.S)<0,-1000000,Social.Welfare.Function))

#Calculate GDP growth rate for both regions
# out$Growth.Rate_N<-c(NA,diff(out$Consumption.N)/out$Consumption.N[1:(length(out$Consumption.N)-1)])
# out$Growth.Rate_S<-c(NA,diff(out$Consumption.S)/out$Consumption.S[1:(length(out$Consumption.S)-1)])
#set Run.ID
# out$Run.ID<-as.numeric(Parameters['Run.ID'])
#set social function value
# out$Social.Welfare.Function<-Social.Welfare.Function
#if (verbose==FALSE) {
#                      return(Social.Welfare.Function)
#                    } else {
                      return(out)
#                           }
 }
