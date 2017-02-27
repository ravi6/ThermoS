within ThermoS.Uops.Valves;
partial model partialValve
import Modelica.Media.IdealGases.Common.MixtureGasNasa.h_TX;

  // A simple Valve
  replaceable package Medium = PartialMixtureMedium;

  // Specify that our Medium is used in/outlet
  FluidPort inlet(redeclare package Medium = Medium);
  FluidPort outlet(redeclare package Medium = Medium);

  // Valve Coefficient   (1 m3/s @ one bar differential with 1kg/m3 density)
  parameter Real cv=1.0/sqrt(1e5);
  parameter Real dpTol = 100 ;  // Pressure Drop Tolerance used for state and flow transition

 // Medium.ThermodynamicState state (T(start=300), p(start=1e5), X(start=Medium.reference_X));
  Medium.BaseProperties med ; // (preferredMediumStates=false) ;  
 
equation

/* No change in Enthalpy in the Valve
    * inStream returns value when flow is into the device
    * Looks trivial for one to one connections. But is designed to handle
    * large number of connection branches without singularity
*/

  outlet.h_outflow  = inStream(inlet.h_outflow); // flow is  inlet to outlet
  inlet.h_outflow   = inStream(outlet.h_outflow); // reverse flow case

  inlet.Xi_outflow  = inStream(outlet.Xi_outflow); // pass composition un_altered
  outlet.Xi_outflow = inStream(inlet.Xi_outflow); // pass composition un_altered

  // Mass balance
  inlet.m_flow + outlet.m_flow = 0;

  med.p  = max(inlet.p, outlet.p);
  med.h  = inlet.h_outflow  ;
  med.Xi = inlet.Xi_outflow ;

end partialValve;
