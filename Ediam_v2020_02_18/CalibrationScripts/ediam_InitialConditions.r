#==========================================================================================================================================================================================================
#Define Initial Conditions
#==========================================================================================================================================================================================================
#subindices are undersores "_", supraindices are periods ".", subindices go before periods

ediamEmpiricalParameters<-function(Parameters,t)
{
#
#Determine initial conditions
 alfa.N<-as.numeric(Parameters['alfa.N'])
 phi.N<-as.numeric(Parameters['phi.N'])
 epsilon.N<-as.numeric(Parameters['epsilon.N'])
 alfa.S<-as.numeric(Parameters['alfa.S'])
 phi.S<-as.numeric(Parameters['phi.S'])
 epsilon.S<-as.numeric(Parameters['epsilon.S'])


#Data
#t<-1
 Y_re0.N <<- Y_re.Nh[t]
 Y_ce0.N <<- Y_ce.Nh[t]
 Y_re0.S <<- Y_re.Sh[t]
 Y_ce0.S <<- Y_ce.Sh[t]
 Price.oil0<<-Price.oil.y[t]
 ReToGDP0.N<<-ReToGDP.Nh[t]
 ReToGDP0.S<<-ReToGDP.Sh[t]
 GDP0.N <<- GDP.Nh[t]
 GDP0.S <<- GDP.Sh[t]
 Re0.N <<- Re.Nh[t]
 Re0.S <<- Re.Sh[t]
 L_0.N <<- L.N.y[t]
 L_0.S <<- L.S.y[t]
 #erase after testing
 pi.N<-1
 pi.S<-1

#estimate all known factors:
#north
  Y_0.N<-(Y_re0.N^((epsilon.N-1)/epsilon.N)+Y_ce0.N^((epsilon.N-1)/epsilon.N))^(epsilon.N/(epsilon.N-1))
  Price_ce.N<-((Y_re0.N/Y_ce0.N)^((epsilon.N-1)/epsilon.N)+1)^(1/(epsilon.N-1))
  Price_re.N<-( (Y_re0.N/Y_ce0.N)^(-1/epsilon.N) )* Price_ce.N

#Empirical parameters-Advanced Region-
 a.alfa_2.N<-1
 b.alfa_2.N<--1*((GDP0.N/Y_0.N)*((Y_ce0.N*(Price_ce.N^epsilon.N))/(Price.oil0*Re0.N ))+2*alfa.N+1)
 c.alfa_2.N<-(1- ( (Price_re.N^(1-epsilon.N)) + (Price_ce.N^(1-epsilon.N)) ) *alfa.N^2)/(Price_ce.N^(1-epsilon.N))
 alfa_2.N<-(-1*b.alfa_2.N - ( (b.alfa_2.N^2)-4*a.alfa_2.N*c.alfa_2.N)^0.5 )/(2*a.alfa_2.N)
#
#second estimate alfa 1
  alfa_1.N<-alfa.N-alfa_2.N
#third estimate price of final good "Y"
  Price_Y0.N<-(1/alfa_2.N) * Price.oil0 * ( Re0.N/Y_ce0.N) * ( (Y_re0.N/Y_ce0.N)^( (epsilon.N-1)/epsilon.N ) + 1 )^(1/(1-epsilon.N))

#============================
#old approach
#first estimate alfa 2
#  a.alfa_2.N<-1
#  b.alfa_2.N<-( (1/ReToGDP0.N) +1 - 2*alfa.N )
#  c.alfa_2.N<-( alfa.N^2 - pi.N )*(
#                                (Y_re0.N/Y_ce0.N)^( (epsilon.N-1)/epsilon.N ) + 1
#                                )

#  alfa_2.N<-(-1*b.alfa_2.N + ( (b.alfa_2.N^2)-4*a.alfa_2.N*c.alfa_2.N)^0.5 )/(2*a.alfa_2.N)
#second estimate alfa 1
#  alfa_1.N<-alfa.N-alfa_2.N
#third estimate price of final good "Y"
#  Price_Y0.N<-(1/alfa_2.N) * Price.oil0 * ( Re0.N/Y_ce0.N) * ( (Y_re0.N/Y_ce0.N)^( (epsilon.N-1)/epsilon.N ) + 1 )^(1/(1-epsilon.N))
#============================


#fourth initial costs of technology and auxiliary
  epsi_re.N<-alfa.N^2
  epsi_ce.N<-alfa_1.N^2
  phi_1.N<-(1-alfa_1.N)*(1-epsilon.N)


#Empirical parameters-Emerging Region-
  Y_0.S<-(Y_re0.S^((epsilon.S-1)/epsilon.S)+Y_ce0.S^((epsilon.S-1)/epsilon.S))^(epsilon.S/(epsilon.S-1))
  Price_ce.S<-( ((Y_re0.S/Y_ce0.S)^((epsilon.S-1)/epsilon.S)+1)^(1/(epsilon.S-1)) )
  Price_re.S<-( (Y_re0.S/Y_ce0.S)^(-1/epsilon.S) )* Price_ce.S

  a.alfa_2.S<-1
  b.alfa_2.S<--1*((GDP0.S/Y_0.S)*((Y_ce0.S*(Price_ce.S^epsilon.S))/(Price.oil0*Re0.S ))+2*alfa.S+1)
  c.alfa_2.S<-(1- ( (Price_re.S^(1-epsilon.S)) + (Price_ce.S^(1-epsilon.S)) ) *alfa.S^2)/(Price_ce.S^(1-epsilon.S))

  alfa_2.S<-(-1*b.alfa_2.S - ( (b.alfa_2.S^2)-4*a.alfa_2.S*c.alfa_2.S)^0.5 )/(2*a.alfa_2.S)

#============================
#old approah
#first estimate alfa 2
#  a.alfa_2.S<-1
#  b.alfa_2.S<-( (1/ReToGDP0.S) +1 - 2*alfa.S )
#  c.alfa_2.S<-( alfa.S^2 - pi.S )*(
#                                (Y_re0.S/Y_ce0.S)^( (epsilon.S-1)/epsilon.S ) + 1
#                                )

#  alfa_2.S<-(-1*b.alfa_2.S + ( (b.alfa_2.S^2)-4*a.alfa_2.S*c.alfa_2.S)^0.5 )/(2*a.alfa_2.S)
#second estimate alfa 1
#  alfa_1.S<-alfa.S-alfa_2.S
#third estimate price of final good "Y"
#  Price_Y0.S<-(1/alfa_2.S) * Price.oil0 * ( Re0.S/Y_ce0.S) * ( (Y_re0.S/Y_ce0.S)^( (epsilon.S-1)/epsilon.S ) + 1 )^(1/(1-epsilon.S))
#=================================

#
 #second estimate alfa 1
   alfa_1.S<-alfa.S-alfa_2.S
 #third estimate price of final good "Y"
   Price_Y0.S<-(1/alfa_2.S) * Price.oil0 * ( Re0.S/Y_ce0.S) * ( (Y_re0.S/Y_ce0.S)^( (epsilon.S-1)/epsilon.S ) + 1 )^(1/(1-epsilon.S))
 #fourth initial costs of technology and auxiliary
   epsi_re.S<-alfa.S^2
   epsi_ce.S<-alfa_1.S^2
   phi_1.S<-(1-alfa_1.S)*(1-epsilon.S)

#return values
EmpiricalParams <- data.frame(
                              alfa_1.N = alfa_1.N,
                              alfa_2.N = alfa_2.N,
                              epsi_re0.N = epsi_re.N,
                              epsi_ce0.N = epsi_ce.N,
                              phi_1.N = phi_1.N,
                              Price.final.N = Price_Y0.N,
                              alfa_1.S = alfa_1.S,
                              alfa_2.S = alfa_2.S,
                              Price.final.S = Price_Y0.S,
                              epsi_re0.S = epsi_re.S,
                              epsi_ce0.S = epsi_ce.S,
                              phi_1.S = phi_1.S,
                              Y_0.N = Y_0.N,
                              Y_0.S =Y_0.S
                            )
return(EmpiricalParams)

}

#subindices are undersores "_", supraindices are periods ".", subindices go before periods
ediamInitialConditions<-function(Parameters)
{
#Determine initial conditions
 alfa.N<-as.numeric(Parameters['alfa.N'])
 phi.N<-as.numeric(Parameters['phi.N'])
 epsilon.N<-as.numeric(Parameters['epsilon.N'])
 alfa.S<-as.numeric(Parameters['alfa.S'])
 phi.S<-as.numeric(Parameters['phi.S'])
 epsilon.S<-as.numeric(Parameters['epsilon.S'])
  Eta_re.N<-as.numeric(Parameters['Eta_re.N'])
 Eta_ce.N<-as.numeric(Parameters['Eta_ce.N'])
 Eta_re.S<-as.numeric(Parameters['Eta_re.S'])
 Eta_ce.S<-as.numeric(Parameters['Eta_ce.S'])
 alfa_1.N <- as.numeric(Parameters['alfa_1.N'])
 alfa_2.N <- as.numeric(Parameters['alfa_2.N'])
 epsi_re.N <- as.numeric(Parameters['epsi_re0.N'])
 epsi_ce.N <- as.numeric(Parameters['epsi_ce0.N'])
 phi_1.N <- as.numeric(Parameters['phi_1.N'])
 Price.final.N <- as.numeric(Parameters['Price.final.N'])
 alfa_1.S <- as.numeric(Parameters['alfa_1.S'])
 alfa_2.S <- as.numeric(Parameters['alfa_2.S'])
 Price.final.S <- as.numeric(Parameters['Price.final.S'])
 epsi_re.S <- as.numeric(Parameters['epsi_re0.S'])
 epsi_ce.S <- as.numeric(Parameters['epsi_ce0.S'])
 phi_1.S <- as.numeric(Parameters['phi_1.S'])
#
#test
Y_0.N <-as.numeric(Parameters['Y_0.N'])
Y_0.S <- as.numeric(Parameters['Y_0.S'])

#Data
Y_re0.N <<- Y_re.Nh[1]
 Y_ce0.N <<- Y_ce.Nh[1]
 Y_re0.S <<- Y_re.Sh[1]
 Y_ce0.S <<- Y_ce.Sh[1]
 Price.oil0<<-Price.oil.y[1]
 ReToGDP0.N<<-ReToGDP.Nh[1]
 ReToGDP0.S<<-ReToGDP.Sh[1]
 GDP0.N <<- GDP.Nh[1]
 GDP0.S <<- GDP.Sh[1]
Re0.N <<- Re.Nh[1]
 Re0.S <<- Re.Sh[1]
 L_0.N <<- L.N.y[1]
 L_0.S <<- L.S.y[1]

 S_0<-as.numeric(Parameters['S_0'])

#
#Initial Productivity conditions are determined by the initial levels of production of energy and initial prices of the exhaustable resource

#fifth productivity carbon sector
  A_ce0.N<-( L_0.N^((1-alfa.N)/(alfa_1.N-1)) ) *
           ( Y_ce0.N^((1-alfa.N)/(1-alfa_1.N)) ) *
           ( (epsi_ce.N/(alfa_1.N^2))^(alfa_1.N/(1-alfa_1.N)) ) *
           ( (Price.oil0/(Price.final.N*alfa_2.N))^(alfa_2.N/(1-alfa_1.N)) ) *
           ( ( (Y_re0.N/Y_ce0.N)^( (epsilon.N-1)/epsilon.N ) + 1 )^( (alfa.N+phi.N)/phi_1.N ) )
#sixth productivity renewable sector
  A_re0.N<- (
              ( ( (epsi_re.N^alfa.N) * (alfa_1.N^(2*alfa_1.N))* (alfa_2.N^alfa_2.N)* (A_ce0.N^(1-alfa_1.N)) )/( (epsi_ce.N^alfa_1.N)*(alfa.N^(2*alfa.N)) ) ) *
             ( ( Price.final.N/Price.oil0 )^alfa_2.N ) *
             ( (Y_re0.N/Y_ce0.N)^ (1/epsilon.N) )
           )^(1/(1-alfa.N))


 #Initial Conditions-Emerging Region-

 #fifth productivity carbon sector
   A_ce0.S<-( L_0.S^((1-alfa.S)/(alfa_1.S-1)) ) *
            ( Y_ce0.S^((1-alfa.S)/(1-alfa_1.S)) ) *
            ( (epsi_ce.S/(alfa_1.S^2))^(alfa_1.S/(1-alfa_1.S)) ) *
            ( (Price.oil0/(Price.final.S*alfa_2.S))^(alfa_2.S/(1-alfa_1.S)) ) *
            ( ( (Y_re0.S/Y_ce0.S)^( (epsilon.S-1)/epsilon.S ) + 1 )^( (alfa.S+phi.S)/phi_1.S ) )
 #sixth productivity renewable sector
   A_re0.S<- (
               ( ( (epsi_re.S^alfa.S) * (alfa_1.S^(2*alfa_1.S))* (alfa_2.S^alfa_2.S)* (A_ce0.S^(1-alfa_1.S)) )/( (epsi_ce.S^alfa_1.S)*(alfa.S^(2*alfa.S)) ) ) *
              ( ( Price.final.S/Price.oil0 )^alfa_2.S ) *
              ( (Y_re0.S/Y_ce0.S)^ (1/epsilon.S) )
             )^(1/(1-alfa.S))

#Determine initial levels of expected profits for each sector:
#Advanced region
  cR0.N<-Price.oil0/Price.final.N
  RelPrice.N<-( (epsi_re.N^alfa.N)*(alfa_1.N^(2*alfa_1.N))*(alfa_2.N^alfa_2.N)*(A_ce0.N^(1-alfa_1.N)) )/
                 ( (cR0.N^alfa_2.N)*(epsi_ce.N^alfa_1.N)*(alfa.N^(2*alfa.N))*(A_re0.N^(1-alfa.N)) )
 RelLabor.N<-((1+0)^epsilon.N)*
            ( (
               ( (alfa.N^(2*alfa.N))*(cR0.N^alfa_2.N)*(epsi_ce.N^alfa_1.N) )/
               ( ((1-0)^alfa.N)*(alfa_2.N^alfa_2.N)*(alfa_1.N^(2*alfa_1.N))*(epsi_re.N^alfa.N) )
               )^(epsilon.N-1) )*
            (
               ( A_re0.N^(-1*phi.N) )/
               ( A_ce0.N^(-1*(1-alfa_1.N)*(1-epsilon.N)) )
             )
#Initial profits clean sector
      Labor_re.N<-(RelLabor.N*L_0.N)/(1+RelLabor.N)
      Price_re.N<-RelPrice.N/(RelPrice.N^(1-epsilon.N)+(1)^(1-epsilon.N))^(1/(1-epsilon.N))
      Xtech_re.N<-((((alfa.N^2)*Price_re.N)/((1)*epsi_re.N))^(1/(1-alfa.N)))*Labor_re.N*A_re0.N
Profits_re0.N<-Eta_re.N*epsi_re.N*((1-alfa.N)/alfa.N)*Xtech_re.N
#Initial profits dirty sector
      Labor_ce.N<-L_0.N/(RelLabor.N+1)
      Price_ce.N<-Price_re.N/RelPrice.N
      Re.N<-(alfa_1.N^2/epsi_ce.N)^( alfa_1.N/(1-alfa.N) )*
          (alfa_2.N*A_ce0.N/cR0.N)^((1-alfa_1.N)/(1-alfa.N))*
          Labor_ce.N*
          Price_ce.N^(1/(1-alfa.N))
      Xtech_ce.N<- (
                               (
                                ( (alfa_1.N^2)*(Price_ce.N)*(Re.N^alfa_2.N) )/
                                ( epsi_ce.N )
                               )^( 1/(1-alfa_1.N) )
                              )*
                              ( Labor_ce.N^( (1-alfa.N)/(1-alfa_1.N) ) )*
                              ( A_ce0.N )
Profits_ce0.N<-Eta_ce.N*epsi_ce.N*((1-alfa_1.N)/alfa_1.N)*Xtech_ce.N
#
#Emerging region
cR0.S<-Price.oil0/Price.final.S
RelPrice.S<-( (epsi_re.S^alfa.S)*(alfa_1.S^(2*alfa_1.S))*(alfa_2.S^alfa_2.S)*(A_ce0.S^(1-alfa_1.S)) )/
               ( (cR0.S^alfa_2.S)*(epsi_ce.S^alfa_1.S)*(alfa.S^(2*alfa.S))*(A_re0.S^(1-alfa.S)) )
RelLabor.S<-((1+0)^epsilon.S)*
          ( (
             ( (alfa.S^(2*alfa.S))*(cR0.S^alfa_2.S)*(epsi_ce.S^alfa_1.S) )/
             ( ((1-0)^alfa.S)*(alfa_2.S^alfa_2.S)*(alfa_1.S^(2*alfa_1.S))*(epsi_re.S^alfa.S) )
             )^(epsilon.S-1) )*
          (
             ( A_re0.S^(-1*phi.S) )/
             ( A_ce0.S^(-1*(1-alfa_1.S)*(1-epsilon.S)) )
           )
#Initial profits clean sector
    Labor_re.S<-(RelLabor.S*L_0.S)/(1+RelLabor.S)
    Price_re.S<-RelPrice.S/(RelPrice.S^(1-epsilon.S)+(1)^(1-epsilon.S))^(1/(1-epsilon.S))
    Xtech_re.S<-((((alfa.S^2)*Price_re.S)/((1)*epsi_re.S))^(1/(1-alfa.S)))*Labor_re.S*A_re0.S
Profits_re0.S<-Eta_re.S*epsi_re.S*((1-alfa.S)/alfa.S)*Xtech_re.S
#Initial profits dirty sector
    Labor_ce.S<-L_0.S/(RelLabor.S+1)
    Price_ce.S<-Price_re.S/RelPrice.S
    Re.S<-(alfa_1.S^2/epsi_ce.S)^( alfa_1.S/(1-alfa.S) )*
        (alfa_2.S*A_ce0.S/cR0.S)^((1-alfa_1.S)/(1-alfa.S))*
        Labor_ce.S*
        Price_ce.S^(1/(1-alfa.S))
    Xtech_ce.S<- (
                             (
                              ( (alfa_1.S^2)*(Price_ce.S)*(Re.S^alfa_2.S) )/
                              ( epsi_ce.S )
                             )^( 1/(1-alfa_1.S) )
                            )*
                            ( Labor_ce.S^( (1-alfa.S)/(1-alfa_1.S) ) )*
                            ( A_ce0.S )
Profits_ce0.S<-Eta_ce.S*epsi_ce.S*((1-alfa_1.S)/alfa_1.S)*Xtech_ce.S

#define consumption
  Consumption.N<-(Y_0.N-epsi_re.N*Xtech_re.N-epsi_ce.N*Xtech_ce.N-cR0.N*Re.N)
  Consumption.S<-(Y_0.S-epsi_re.S*Xtech_re.S-epsi_ce.S*Xtech_ce.S-cR0.S*Re.S)
  GDP.N<-Price.final.N * Consumption.N
  GDP.S<-Price.final.S * Consumption.S



#Define vector of intital conditions
  InitialConditions <-list(
                        StateVariables = c( A_re = max(A_re0.N,A_re0.S),
                                            A_ce = max(A_ce0.N,A_ce0.S),
                                            A_re.N = A_re0.N,
                                            A_ce.N = A_ce0.N,
                                            A_re.S = A_re0.S,
                                            A_ce.S = A_ce0.S,
                                            E_re.N = Xtech_re.N,
                                            E_ce.N = Xtech_ce.N,
                                            E_re.S = Xtech_re.S,
                                            E_ce.S = Xtech_ce.S,
                                            S= S_0),
                        Parameters = c(   Profits0.N = min(Profits_re0.N,Profits_ce0.N),
                                          Profits0.S = min(Profits_re0.S,Profits_ce0.S),
                                          Xtech_re0.N = Xtech_re.N,
                                          Xtech_ce0.N = Xtech_ce.N,
                                          Xtech_re0.S= Xtech_re.S,
                                          Xtech_ce0.S = Xtech_ce.S,
                                          GDP.N=GDP.N,
                                          GDP.S=GDP.S)
                                      )
  return(InitialConditions)
}
