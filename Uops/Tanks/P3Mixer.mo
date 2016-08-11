// Version 1.0
// Date:   9th December 2013
within ThermoS.Uops.Tanks;
model P3Mixer
/*
*  A Three Port Mixer ... for testing
*/

    replaceable package Medium = PartialMixtureMedium ;
    FluidPort port1(redeclare package Medium = Medium) ;
    FluidPort port2(redeclare package Medium = Medium) ;
    FluidPort port3(redeclare package Medium = Medium) ;

//  Parameters
  parameter    Volume   	   vol   = 1    ;   // Tank Volume (m3)

// State Variables
    Mass			m		;	// Mass of Gas in the vessal 
    Medium.Temperature		T		;
    Medium.AbsolutePressure	p		;
    Medium.MassFraction		Xi[Medium.nXi]	;
    Medium.BaseProperties       medium          ;

  equation

    medium.T = T ;
    medium.p = p ;
    medium.Xi = Xi ;


     m = medium.d * vol ; 

     // Mass and Component Balance
     der(m) = port1.m_flow + port2.m_flow + port3.m_flow   ;  // Mass Balance

     // Component Balance
     der(Xi*m) = actualStream(port1.Xi_outflow) * port1.m_flow 
              +  actualStream(port2.Xi_outflow) * port2.m_flow 
              +  actualStream(port3.Xi_outflow) * port3.m_flow ;

     // Enthalpy Balance
     der(m*medium.h) =  port1.m_flow * actualStream(port1.h_outflow)
                     +  port2.m_flow * actualStream(port2.h_outflow)
                     +  port3.m_flow * actualStream(port3.h_outflow)
    	             + vol * der(p) ; 

     // Assume gas in tank is well mixed (ie. its contents are at outlet condition)
        port1.Xi_outflow = medium.Xi ; port2.Xi_outflow = medium.Xi ; port3.Xi_outflow = medium.Xi ;
	port1.h_outflow  = medium.h  ; port2.h_outflow  = medium.h  ; port3.h_outflow  = medium.h  ;
        port1.p = p ; port2.p = p ; port3.p = p ;

end P3Mixer;
