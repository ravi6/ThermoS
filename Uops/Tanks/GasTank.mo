// Version 1.0
// Date:   9th December 2013
within ThermoS.Uops.Tanks;
model GasTank
/*
*  A Gas Storage Vessel with two ports
*     Note inlet and outlet
*      flow can eventuate in any direction or nominal names
*/

    replaceable package Medium = PartialMixtureMedium ;
    FluidPort inlet(redeclare package Medium = Medium) ;
    FluidPort outlet(redeclare package Medium = Medium) ;  

//  Parameters
  parameter    Volume   	   vol   = 10    ;   // Tank Volume (m3)
  parameter    EnthalpyFlowRate    Q_in  = 0    ;

// State Variables
    Mass			m		;	// Mass of Gas in the vessal 
//    Medium.ThermodynamicState	state		;		
    Medium.Temperature		T		;
    Medium.AbsolutePressure	p		;
    Medium.MassFraction		Xi[Medium.nXi]	;
    Medium.SpecificEnthalpy	h		;
    
    Medium.BaseProperties  medium ; //(p(start=1e5), T(start=300),
                                   // Xi(start={1/Medium.nXi, 1/Medium.nXi})) ;

  equation
 //     state = Medium.setState_pTX(p, T, Xi) ;
 //     h     = Medium.specificEnthalpy(state) ;

    medium.T = T ; medium.p = p ; medium.Xi = Xi ;
    medium.h = h ;
     

//     m = Medium.density(state) * vol ; // Mass and Component Balance
     m = medium.d * vol ; // Mass and Component Balance
     der(m) = inlet.m_flow + outlet.m_flow   ;  // Mass Balance

     // Component Balance
     der(Xi*m) = actualStream(inlet.Xi_outflow) * inlet.m_flow 
                  + actualStream(outlet.Xi_outflow) * outlet.m_flow ;

     // Enthalpy Balance
     der(m*h) = Q_in + inlet.m_flow * actualStream(inlet.h_outflow)
                     + outlet.m_flow * actualStream(outlet.h_outflow)
		     + vol * der(p) ; 

     // Assume gas in tank is well mixed (ie. its contents are at outlet condition)
        // state = inlet.Medium.setState_phX(p, h, Xi) ;  // gives 5 equations
        inlet.Xi_outflow = Xi ;
	inlet.h_outflow  = h  ;
	outlet.Xi_outflow = Xi ;
        outlet.h_outflow = h ;

     // No pressure drop acroos the unit
      inlet.p = outlet.p ;
      inlet.p = p ;
end GasTank;
