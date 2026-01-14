// Version 1.0
// Date:   9th December 2013
within ThermoS.Uops.Tanks;
model GasTank
/*
*  A Gas Storage Vessel with two ports
*     Note inlet and outlet
*      flow can eventuate in any direction 
*      These model eqns. ensure static State selection
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

    Medium.BaseProperties       medium        ;

    Mass                 M               ; // Gas mass in tank
    InternalEnergy       U               ; // of tank contents
    
  equation

    // Tank p,T for ease of access
    T = medium.T  ;
    p = medium.p  ;
    Xi = medium.Xi ;
    outlet.Xi_outflow = medium.Xi ;
    
    // Mass Balance
    M = medium.d * vol ;
    der(M) = inlet.m_flow + outlet.m_flow ;

    // Component Balance
    der (M * Xi) = actualStream(inlet.Xi_outflow) * inlet.m_flow 
                      + actualStream(outlet.Xi_outflow) * outlet.m_flow ;

     // Enthalpy Balance
     U = M * medium.u ;
     vol * der(U)  = - Q_loss 
                     + inlet.m_flow * actualStream(inlet.h_outflow)
                     + outlet.m_flow * actualStream(outlet.h_outflow);

     // Assume gas in tank is well mixed 
         inlet.Xi_outflow = outlet.Xi_outflow ;
	 inlet.h_outflow  = outlet.h_outflow  ;
         outlet.h_outflow = medium.h ;

     // No pressure drop acroos the unit
      inlet.p = outlet.p ;
      inlet.p = medium.p ;
end GasTank;
