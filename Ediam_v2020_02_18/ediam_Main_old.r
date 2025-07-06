## =======================================================================================================================================================================================================
## Exploratory Dynamic Integreated Assessment Model
## =======================================================================================================================================================================================================
ediam<-function (policies, dir.model, verbose=FALSE)
{
#==========================================================================================================================================================================================================
#Define policy parameters
#==========================================================================================================================================================================================================
 policy.vector<-c(
    #Carbon tax
	    ce.tax.N=signif(as.numeric(policies[1]), digits = 2),
      ce.tax.S=signif(as.numeric(policies[2]), digits = 2),
    #Technology push in North
	    Tec.subsidy.N = signif(as.numeric(policies[3]), digits = 2),
	    RD.subsidy.N= signif(as.numeric(policies[4]), digits = 2),
    #Traditional Green Climate Fund
	    Tec.subsidy.S = signif(as.numeric(policies[5]), digits = 2),
 	    Tec.subsidy.GF.N = signif(as.numeric(policies[6]), digits = 2),
	  #R&D Green Climate Fund
	    RD.subsidy.S = signif(as.numeric(policies[7]), digits = 2),
	    RD.subsidy.GF.N =signif(as.numeric(policies[8]), digits = 2),
      policy.half.life = as.numeric(policies[9]))
#==========================================================================================================================================================================================================

#==========================================================================================================================================================================================================
#Define model parameters
#==========================================================================================================================================================================================================
Parameters <- c(Run.ID = 1, #as.numeric(params['Run.ID'])
                EndTime = 300, #as.numeric(params['EndTime'])
                TimeStep = 5, #as.numeric(params['TimeStep'])
                alfa = 0.33, #as.numeric(params['alfa']),
                alfa_1= 0.5*0.33,
                alfa_2= 0.5*0.33,
                epsilon = 5.0,#as.numeric(params['epsilon']),
                Gamma_re = 0.25 ,#as.numeric(params['Gamma_re']),
                Gamma_ce = 0.25,#as.numeric(params['Gamma_ce']),
                Eta_re= 0.02,#as.numeric(params['Eta_re']),
                Eta_ce= 0.02,#as.numeric(params['Eta_ce']),
                Nu_re = 0.02,#as.numeric(params['Nu_re']),
                Nu_ce= 0.02,#as.numeric(params['Nu_ce']),
                qsi= as.numeric(0.0100539),#as.numeric(params['qsi']),
                Delta.S = 0.001822767,#as.numeric(params['Delta.S']),
				        Delta.Temp.Disaster = as.numeric(6.0) ,#as.numeric(params['Delta.Temp.Disaster']),
				        Beta.Delta.Temp = 5.0 ,#as.numeric(params['Beta.Delta.Temp']),
				        CO2.base = 289.415,#as.numeric(params['CO2.base']),
				        labor.growth.N = 0.0 ,#as.numeric(params['labor.growth_N']),
				        labor.growth.S = 0.0 ,#as.numeric(params['labor.growth_S']),
				        rho = 0.008,#as.numeric(params['rho']),
				        lambda.S = 0.1443,#as.numeric(params['lambda.S']),
				        sigma.utility = 2.0,#as.numeric(params['sigma.utility']),
                Y_re0.N = 45.55, #as.numeric(params['Yre.0_N']),
                Y_ce0.N = 193.2, #as.numeric(params['Yce.0_N'])
                Y_re0.S = 27.82, #as.numeric(params['Yre.0_S'])
                Y_ce0.S = 257.54, #as.numeric(params['Yce.0_S'])
                CO2.Concentration_0 = 382.2461 #as.numeric(params['CO2.Concentration_0'])
                )
#Add additional parameters
  Parameters<-c(Parameters,
                phi = (1-as.numeric(Parameters['alfa']))*(1-as.numeric(Parameters['epsilon'])),
                phi_1 = (1-as.numeric(Parameters['alfa_1']))*(1-as.numeric(Parameters['epsilon'])),
                epsi = as.numeric(Parameters['alfa'])^2,
                cR_0 = 1.0,
                L_0.N = 1.0,
                L_0.S = 4.0,
                policy.vector)

#==========================================================================================================================================================================================================

#==========================================================================================================================================================================================================
#Define Initial Conditions
#==========================================================================================================================================================================================================

#Load Initial Conditions Function
 source(paste(dir.model,"ediam_InitialConditions.r",sep=""))

#Determine initial conditions
 InitialConditions <-ediamInitialConditions(as.numeric(Parameters['alfa']),
                                             as.numeric(Parameters['alfa_1']),
                                             as.numeric(Parameters['alfa_2']),
                                             as.numeric(Parameters['phi']),
                                             as.numeric(Parameters['phi_1']),
                                             as.numeric(Parameters['epsi']),
                                             as.numeric(Parameters['epsilon']),
                                             as.numeric(Parameters['L_0.N']),
                                             as.numeric(Parameters['L_0.S']),
                                             as.numeric(Parameters['Y_re0.N']),
                                             as.numeric(Parameters['Y_ce0.N']),
                                             as.numeric(Parameters['Y_re0.S']),
                                             as.numeric(Parameters['Y_ce0.S']),
                                             as.numeric(Parameters['cR_0']),
                                             as.numeric(Parameters['CO2.Concentration_0']))

#==========================================================================================================================================================================================================
#Load EDIAM Differential Equations Structure and Execute Simulation
#==========================================================================================================================================================================================================

#Load EDIAM Differential Equations
 source(paste(dir.model,"ediam_Equations.r",sep=""))

#Load Data Series
#for Oil
 cR.Data <<- approxfun( x = seq(0,300,5), #times
                       y = rep(1,61), # Actual data points
                       method = "linear",
                       rule = 2)
#for Population

#Run the model
 out <- as.data.frame(
                      ode(y = InitialConditions,
                          times = seq(0,as.numeric(Parameters['EndTime']),as.numeric(Parameters['TimeStep'])),
                          func = EdiamEquations,
                          parms = Parameters,
                          method ="rk4" ))

#==========================================================================================================================================================================================================
# Post-processing data and validation of constraints
#==========================================================================================================================================================================================================

#Define objetive function
 Social.Welfare.Function<-sum(as.numeric(out$Utility.Consumer_N))+sum(as.numeric(out$Utility.Consumer_S))

#Check constraints are valid
#Consumption constraint -In no year consumption can be negative-
 Social.Welfare.Function<-ifelse(min(out$Consumption.N)<0,-1000000,
                                 ifelse(min(out$Consumption.S)<0,-1000000,Social.Welfare.Function))

#Budget constraint -the cost of policies cannot exceed total budget-
 Social.Welfare.Function<-ifelse(sum(out$Budget.function.N)<0,-1000000,
                                 ifelse(sum(out$Budget.function.S)<0,-1000000,Social.Welfare.Function))

#Calculate consumption growth rate for both regions
 out$Growth.Rate_N<-c(NA,diff(out$Consumption.N)/out$Consumption.N[1:(length(out$Consumption.N)-1)])
 out$Growth.Rate_S<-c(NA,diff(out$Consumption.S)/out$Consumption.S[1:(length(out$Consumption.S)-1)])
#set Run.ID
 out$Run.ID<-as.numeric(Parameters['Run.ID'])
#set social function value
 out$Social.Welfare.Function<-Social.Welfare.Function
if (verbose==FALSE) {
                      return(Social.Welfare.Function)
                    } else {
                      return(out)
                           }
 }
