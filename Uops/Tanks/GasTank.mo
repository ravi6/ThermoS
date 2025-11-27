// Version 1.0
// Date:   9th December 2013
within ThermoS.Uops.Tanks;
model GasTank
/*
*  A Gas Storage Vessel with two ports
*     Note inlet and outlet
*      flow can eventuate in any direction or nominal names
*     Also see how derivatives contain states related medium
*     thus avoiding intermediate variables like, m, h, Xi
*     that can introduce algebraic loops
*/

    replaceable package Medium = PartialMixtureMedium ;
    FluidPort inlet(redeclare package Medium = Medium) ;
    FluidPort outlet(redeclare package Medium = Medium) ;  

//  Parameters
  parameter    Volume   	vol = 1    ;   // Tank Volume (m3)

// State Variables
    Medium.Temperature		T		;
    Medium.AbsolutePressure	p		;
    Medium.MassFraction         Xi[Medium.nXi]  ;
    EnthalpyFlowRate            Q_loss          ; // Heatloss from tank

    Medium.BaseProperties       medium (preferredMediumStates = true) ;
                                   
  equation

    // Tank p,T for ease of access
    medium.T = T ;
    medium.p = p ;
    medium.Xi = outlet.Xi_outflow ;
    Xi = outlet.Xi_outflow ;
 
    // Mass Balance
    vol * der(medium.d) = inlet.m_flow + outlet.m_flow ;

    // Component Balance
     vol * der(medium.d * medium.Xi) = actualStream(inlet.Xi_outflow) * inlet.m_flow 
                  + actualStream(outlet.Xi_outflow) * outlet.m_flow ;

     // Enthalpy Balance
     vol * der(medium.d * medium.h)  = - Q_loss + inlet.m_flow * actualStream(inlet.h_outflow)
                     + outlet.m_flow * actualStream(outlet.h_outflow)
		     + vol * der(medium.p) ; 

     // Assume gas in tank is well mixed 
         inlet.Xi_outflow = outlet.Xi_outflow ;
	 inlet.h_outflow  = outlet.h_outflow  ;
         outlet.h_outflow = medium.h ;

     // No pressure drop acroos the unit
      inlet.p = outlet.p ;
      inlet.p = medium.p ;
end GasTank;
