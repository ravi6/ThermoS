// Version 1.0
// Date:   9th December 2013
within ThermoS.Uops.Tanks;
model Cvolume
/*
*  An adiabatic Volume with one port 
*/

    replaceable package Medium = PartialMixtureMedium ;
    FluidPort port(redeclare package Medium = Medium) ;

//  Parameters
  parameter    Volume   	   vol   = 1    ;   // Tank Volume (m3)

// State Variables
    Mass			m		;	// Mass of Gas in the vessal 
    Medium.ThermodynamicState	state		;		
    Medium.Temperature		T		;
    Medium.AbsolutePressure	p		;
    Medium.MassFraction		Xi[Medium.nXi]	;
    Medium.SpecificEnthalpy	h		;

  equation
      state = Medium.setState_pTX(p, T, Xi) ;
      h     = Medium.specificEnthalpy(state) ;

     m = Medium.density(state) * vol ; 

     // Mass and Component Balance
     der(m) = port.m_flow   ;  // Mass Balance

     // Component Balance
     der(Xi*m) = actualStream(port.Xi_outflow) * port.m_flow ;

     // Enthalpy Balance
     der(m*h) =  port.m_flow * actualStream(port.h_outflow)
		     + vol * der(p) ; 

     // Assume gas in tank is well mixed (ie. its contents are at outlet condition)
        port.Xi_outflow = Xi ;
	port.h_outflow  = h  ;
        port.p = p ;

end Cvolume;
