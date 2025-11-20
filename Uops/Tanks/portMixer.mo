// Version 1.0
// Date:   9th December 2013
within ThermoS.Uops.Tanks;
model portMixer
/*
*  A Three Port Mixer  with finite volume
*/

    replaceable package Medium = PartialMixtureMedium ;

//  Parameters
  parameter    Volume   	   vol   	= 1e-6 ;   // Tank Volume (m3)
  parameter    Boolean             Adiabatic  	= false ;   // default is Isothermal
  parameter    Medium.Temperature  Tset 	= 300   ;   //  at 300K
  parameter    Integer		   N 		= 2     ;   // Number of Ports

  FluidPort port[N] (redeclare each package Medium = Medium) ;

// State Variables
    Mass			m (start=1e-6, stateSelect = StateSelect.avoid)  ;// Mass of Gas in the vessal 
    Medium.Temperature		T (start=300)   ;
    Medium.AbsolutePressure	p (start=1e5, stateSelect = StateSelect.prefer)   ;
    Medium.BaseProperties       medium          ;
    HeatFlowRate		Q_in		;

  equation

    if (Adiabatic) then
      Q_in = 0 ;
    else
      T = Tset ;
    end if;

    medium.T = T  ;
    medium.p = p   ;  
    m = medium.d * vol ; 

     // Mass and Component Balance
     der(m) = sum(port[:].m_flow)   ;  // Mass Balance

     // Component Balance
     for j in 1:Medium.nXi loop
       der(m * medium.Xi[j]) = 
              sumi (actualStream(port[i].Xi_outflow[j]) * port[i].m_flow
                                   for i in 1:N)  ;
     end for;

     // Enthalpy Balance
     der(m * medium.h) =  sum (port[i].m_flow * actualStream(port[i].h_outflow)
                                  for i in 1:N)
    	               + vol * der(p)  + Q_in; 

     // Assume gas in tank is well mixed (ie. its contents are at outlet condition)
      for i in 1:N loop
        port[i].Xi_outflow = medium.Xi ; 
	port[i].h_outflow  = medium.h  ;
        port[i].p = p ;
      end for;

end portMixer;
