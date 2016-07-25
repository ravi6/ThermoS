model plant
/*
  Author: Ravi Saripalli
*/


  import Modelica.SIunits.*;
  package Medium  = Modelica.Media.IdealGases.MixtureGases.SimpleNaturalGas;

  Medium.ThermodynamicState state ;
  Temperature T ;
  SpecificEnthalpy h(start=1e5);

  
equation
// Given Enthalpy find T of a Mix
     der(h) = 2.0e5 ;
     state.T = T ;
     //state = Medium.setState_phX(1e5, h, Medium.reference_X[1:2]);
     state = Medium.setState_phX(1e5, h, {0.1, 0.1, 0.1, 0.1, 0.1});
end plant;
