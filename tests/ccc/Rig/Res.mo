within Rig;

model Res
  // A large reservoir   (p, T)
  // Permits reverse flow through its port

  replaceable package Medium = PartialMixtureMedium ;

  // Specify that our Medium is used in outlet
  FluidPort port (redeclare package Medium = Medium)  ; 

  parameter Medium.AbsolutePressure     p               ;
  parameter Medium.Temperature          T               ;

  Medium.BaseProperties  medium  ;

  equation

    medium.T = T ;
    medium.p = p ;
    medium.Xi = Medium.reference_X [1:Medium.nXi] ;
    port.h_outflow = medium.h ;
    port.p = p ;
    port.Xi_outflow = medium.Xi ;

end Res;
