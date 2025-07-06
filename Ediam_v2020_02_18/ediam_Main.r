## =======================================================================================================================================================================================================
## Exploratory Dynamic Integreated Assessment Model
## =======================================================================================================================================================================================================
ediamMain<-function (calib.params, verbose=FALSE)
{
#==========================================================================================================================================================================================================
#Define policy parameters
#==========================================================================================================================================================================================================
 policy.vector<-c(
    #Advanced Region
      ce.tax.N = max(as.numeric(policy.params['ce.tax.N']),0),
       Schedule.ce.tax.N = as.numeric(policy.params['Schedule.ce.tax.N']),
      Tec.subsidy.N = max(min(as.numeric(policy.params['Tec.subsidy.N']),0.90),0),
       Schedule.Tec.subsidy.N = as.numeric(policy.params['Schedule.Tec.subsidy.N']),
      RD.subsidy.N = max(as.numeric(policy.params['RD.subsidy.N']),0),
       Schedule.RD.subsidy.N = as.numeric(policy.params['Schedule.RD.subsidy.N']),
       #policy.half.life.N = as.numeric(policy.params['policy.half.life.N']),
    #Emerging Region
      ce.tax.S= max(as.numeric(policy.params['ce.tax.S']),0),
       Schedule.ce.tax.S = as.numeric(policy.params['Schedule.ce.tax.S']),
      Tec.subsidy.S = max(min(as.numeric(policy.params['Tec.subsidy.S']),0.90),0),
       Schedule.Tec.subsidy.S = as.numeric(policy.params['Schedule.Tec.subsidy.S']),
      RD.subsidy.S = max(as.numeric(policy.params['RD.subsidy.S']),0),
       Schedule.RD.subsidy.S = as.numeric(policy.params['Schedule.RD.subsidy.S'])
      #policy.half.life.S = as.numeric(policy.params['policy.half.life.S'])
                 )
#==========================================================================================================================================================================================================

#==========================================================================================================================================================================================================
#Define model parameters
#==========================================================================================================================================================================================================
 Parameters <- c(Run.ID = 1,
                 EndTime = 87, #simulation runs from 2014 to 2100
                 TimeStep = TimeStep,
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
                 Gamma_re = as.numeric(calib.params['Gamma_re'])*TimeStep,
                 Gamma_ce = as.numeric(calib.params['Gamma_ce'])*TimeStep,
                 Eta_re.N= as.numeric(calib.params['Eta_re.N']),
                 Eta_ce.N= as.numeric(calib.params['Eta_ce.N']),
                 Eta_re.S = as.numeric(calib.params['Eta_re.S']),
                 Eta_ce.S= as.numeric(calib.params['Eta_ce.S']),

                #Environmental parameters
                 earth.atm.weigth = 5.3e18,#[kg]
                 molar.mass.CO2 = 44.0095, #[g/mol]
                 molar.mass.air = 28.97, # [g/mol]
                 qsi= as.numeric(climate.params['qsi']),  # CO2 in ppm per 1 MTOE of FOSSILTPES
                 Delta.S = as.numeric(climate.params['Delta.S'])*TimeStep,# atmophere's absortive capacity
				         Delta.Temp.Disaster = as.numeric(climate.params['Delta.Temp.Disaster']), #temperature at which environmental disater is reached
				         Beta.Delta.Temp = as.numeric(climate.params['Beta.Delta.Temp']),# effect of increase concentrations of co2 on temperature rise
				         CO2.base = as.numeric(climate.params['CO2.base']), # co2 concentrations in pre-industrial times
                 CO2.Disaster = as.numeric(climate.params['CO2.Disaster']), #co2 concentration at which of disaster threshold
                 DeltaTempConstant = as.numeric(climate.params['DeltaTempConstant']), #constant of linear model of temperature rise
				         rho = 0.001,
				         lambda.S = 0.1443*4.0,
				         sigma.utility = 2.0,
                 S_0 = as.numeric(climate.params['S_0'])    # co2 concentration level at initial state of simulation, this is "2014" S_0=CO2.Disaster-CO2[2014]
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
# Input empirical parameters estimated in calibration
# EmpiricalParams<-c(alfa_1.N = 0.29965900 ,
#                    alfa_2.N = 0.02234100 ,
#                    epsi_re0.N = 0.10368400,
#                    epsi_ce0.N  = 0.08981396,
#                    phi_1.N = -0.90203921,
#                    Price.final.N  = 6.51770721,
#                    alfa_1.S = 0.25027601,
#                    alfa_2.S = 0.07272399,
#                    Price.final.S = 2.09667144,
#                    epsi_re0.S = 0.10432900,
#                    epsi_ce0.S = 0.06283066,
#                    phi_1.S = -1.50844467
#                  )

EmpiricalParams<-ediamEmpiricalParameters(Parameters,1)
EmpiricalParams<-apply(EmpiricalParams,2,mean)
#Add to parameters
  Parameters<-c(Parameters,EmpiricalParams)
#Use estimated parameters for setting initial conditions
  InitialConditions <-ediamInitialConditions(Parameters)

#Add additional parameters
  Parameters<-c(Parameters,InitialConditions$Parameters)

#Initial conditions for state variables
  InitialConditions <- InitialConditions$StateVariables

#==========================================================================================================================================================================================================
#Load EDIAM Difference Equations Structure and Execute Simulation
#==========================================================================================================================================================================================================

#Load Data Series
#for Oil
 Price.oil.Data <<- approxfun( x = oil.times, #times
                       y = Price.oil.y, # Actual data points
                       method = "linear",
                       rule = 2)

#for Population
 L.N.Data<<-approxfun( x = pop.times, #times
                       y = L.N.y, # Actual data points
                       method = "linear",
                       rule = 2)
#
 L.S.Data<<-approxfun( x = pop.times, #times
                      y = L.S.y, # Actual data points
                      method = "linear",
                      rule = 2)
#Run the model
 out <- Ediam(end.time = Parameters['EndTime'] ,
              time.step = Parameters['TimeStep'] ,
              y = InitialConditions,
              func = EdiamEquations,
              Parameters = Parameters)

#==========================================================================================================================================================================================================
# Post-processing data and validation of constraints
#==========================================================================================================================================================================================================

#Define objetive function
 Social.Welfare.Function<-sum(as.numeric(out$D.Utility.Consumer.N))+sum(as.numeric(out$D.Utility.Consumer.S))

#Check constraints are valid
#Consumption constraint -In no year consumption can be negative-
 Social.Welfare.Function<-ifelse(min(out$Consumption.N)<0,-100,
                                 ifelse(min(out$Consumption.S)<0,-100,Social.Welfare.Function))

#Budget constraint -the cost of policies cannot exceed total budget-
 Social.Welfare.Function<-ifelse(sum(out$Budget.function.N)<0,-100,
                                 ifelse(sum(out$Budget.function.S)<0,-100,Social.Welfare.Function))

#Calculate GDP growth rate for both regions
 out$GrowthRate.N<-c(NA,diff(out$GDP.N)/out$GDP.N[1:(length(out$GDP.N)-1)])
 out$GrowthRate.S<-c(NA,diff(out$GDP.S)/out$GDP.S[1:(length(out$GDP.S)-1)])


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
