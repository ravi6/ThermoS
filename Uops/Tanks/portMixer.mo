// Version 1.0
// Date:   9th December 2013
within ThermoS.Uops.Tanks;
model portMixer
/*
*  A Three Port Mixer  with finite volume
*/

    replaceable package Medium = PartialMixtureMedium ;

//  Parameters
  parameter    Volume   	   vol   	= 1     ;   // Tank Volume (m3)
  parameter    Boolean             Adiabatic  	= false ;   // default is Isothermal
  parameter    Medium.Temperature  Tset 	= 300   ;   //  at 300K
  parameter    Integer		   N 		= 2     ;   // Number of Ports

  FluidPort port[N] (redeclare each package Medium = Medium) ;

// State Variables
    Mass			m		;	// Mass of Gas in the vessal 
    Medium.Temperature		T		;
    Medium.AbsolutePressure	p		;
    Medium.MassFraction		Xi[Medium.nXi]	;
    Medium.BaseProperties       medium          ;
    Heat			Q_in		;

  equation

    if (Adiabatic) then
      Q_in = 0 ;
    else
      T = Tset ;
    end if;

    medium.T = T ;
    medium.p = p ;
    medium.Xi = Xi ;

     m = medium.d * vol ; 

     // Mass and Component Balance
     der(m) = sum(port[:].m_flow)   ;  // Mass Balance

     // Component Balance
     for j in 1:Medium.nXi loop
       der(Xi[j]*m) = sum(actualStream(port[i].Xi_outflow[j]) * port[i].m_flow
                           for i in 1:N)  ;
     end for;

     // Enthalpy Balance
     der(m*medium.h) =  sum(port[i].m_flow * actualStream(port[i].h_outflow)
                             for i in 1:N)
    	               + vol * der(p)  + Q_in; 

     // Assume gas in tank is well mixed (ie. its contents are at outlet condition)
      for i in 1:N loop
        port[i].Xi_outflow = medium.Xi ; 
	port[i].h_outflow  = medium.h  ;
        port[i].p = p ;
      end for;

end portMixer;
