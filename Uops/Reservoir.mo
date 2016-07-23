within ThermoS.Uops;
model Reservoir
  // A large Reservoir with constant (p, T, & Xi)
  // Permits reverse flow through its port

  replaceable package Medium = PartialMixtureMedium ;
  FluidPort port (redeclare package Medium = Medium)  ; 
	// Specify that our Medium is used in outlet

  parameter Medium.AbsolutePressure     p               ;
  parameter Medium.Temperature          T               ;
  parameter Medium.MassFraction         Xi[Medium.nXi]  ;

  Medium.BaseProperties       medium          ;

//  Medium.ThermodynamicState state ;
//  Medium.ThermodynamicState state (T(start=300), X(start=Medium.reference_X), p(start=1e5))   ;  // This is how you can have implicit start

  equation
//    state 		= Medium.setState_pTX( p, T, Xi ); 
//    port.h_outflow 	= Medium.specificEnthalpy(state) ;

    medium.T = T ;
    medium.p = p ;
    medium.Xi = Xi ;

    port.h_outflow = medium.h ;
    port.Xi_outflow 	= Xi;
    port.p		= p ;

end Reservoir;
