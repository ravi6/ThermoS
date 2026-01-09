// Version 1.0
// Date:   9th December 2013
within ThermoS.Uops.Tanks;
model Cvolume
/*
*  An adiabatic/Isothermal Volume with one port 
*/

    replaceable package Medium = PartialMixtureMedium ;
    FluidPort port(redeclare package Medium = Medium) ;

//  Parameters
  parameter    Volume   	   vol   = 1    ;   // Tank Volume (m3)
  parameter    Boolean   Adiabatic = false      ;   // default Isothermal
  parameter    Medium.Temperature  Tset = 300   ;   // Isothermal Temperature

    Mass			M		;	// Mass of Gas in the vessal 
    Mass			Mx[Medium.nXi]	;	// Mass of each Component
    Enthalpy                    H               ;       // Enthalpy of contents

    Medium.ThermodynamicState	state		;		
    Medium.Temperature		T		;
    Medium.AbsolutePressure	p		;
    Medium.MassFraction		Xi[Medium.nXi]	;
    Medium.SpecificEnthalpy	h		;
    Heat			Q_in            ;

  equation

      if Adiabatic then
         Q_in = 0.0 ;
      else
         T = Tset ;
      end if;

      state = Medium.setState_pTX(p, T, Xi) ;
      h     = Medium.specificEnthalpy(state) ;

     M = Medium.density(state) * vol ; 
     Mx = M * Xi ;

     // Mass and Component Balance
     der(M) = port.m_flow   ;  // Mass Balance

     // Component Balance
     
     der(Mx) = actualStream(port.Xi_outflow) * port.m_flow ;

     // Enthalpy Balance
     H = m * h ;
     der(H) =  port.m_flow * actualStream(port.h_outflow)
		     + vol * der(p)  + Q_in; 

     // Assume gas in tank is well mixed (ie. its contents are at outlet condition)
        port.Xi_outflow = Xi ;
	port.h_outflow  = h  ;
        port.p = p ;

end Cvolume;
