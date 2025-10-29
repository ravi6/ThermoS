within ThermoS.Uops;
model Reservoir
  // A large Reservoir with constant (p, T, & Xi)
  // Permits reverse flow through its port

  replaceable package Medium = PartialMixtureMedium ;
  FluidPort port (redeclare package Medium = Medium)  ; 
	// Specify that our Medium is used in outlet

  parameter Medium.AbsolutePressure     p               ;
  parameter Medium.Temperature          T               ;
  parameter Medium.MassFraction  Xi[Medium.nXi]         ;

  Medium.BaseProperties  medium  ;

  equation

    medium.T = T ;
    medium.p = p ;
    medium.Xi = Xi ;
    port.h_outflow = medium.h ;
    port.p = p ;
    port.Xi_outflow = Xi ;
    port.h_outflow = medium.h ;
   
/*
    port.Xi_outflow = noEvent (if (port.m_flow < 0) then Xi
                      else inStream (port.Xi_outflow)) ;
    port.h_outflow = noEvent (if (port.m_flow < 0) then medium.h 
                      else inStream (port.h_outflow)) ;
*/
end Reservoir;
