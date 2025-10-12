within ThermoS.Uops;
model Feed

  // A material stream feeder
  replaceable package Medium = PartialMixtureMedium ;
  FluidPort outlet (redeclare package Medium = Medium)  ; 

  // Specify that our Medium is used in outlet
  Medium.MassFlowRate    mdot(min=0) 	; 
  Medium.Temperature     T              ;
  Medium.MassFraction    Xi[Medium.nXi]  ;
  Medium.BaseProperties  medium ;

  equation
    medium.p = outlet.p ;
    medium.T =  T ;
    medium.Xi = Xi ;

    outlet.h_outflow = medium.h ;
    outlet.m_flow = - mdot ; // It is negative because it is going out 
    outlet.Xi_outflow = Xi;  
end Feed;
