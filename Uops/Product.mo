within ThermoS.Uops;
model Product
  replaceable package Medium = PartialMixtureMedium;

  FluidPort inlet(redeclare package Medium = Medium);
  // Specify that our Medium is used in outlet

  Medium.AbsolutePressure p;
  Medium.Temperature T;
  Medium.MassFraction Xi[Medium.nXi];
  Medium.BaseProperties  medium   ;

equation
    medium.p = inlet.p ;
    medium.T = T ;
    medium.Xi = Xi ;

  inlet.h_outflow = medium.h ;  
  inlet.Xi_outflow = Xi;  // state X size is alwasy nS

end Product;
