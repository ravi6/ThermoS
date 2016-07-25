model plant
/*
  Author: Ravi Saripalli
*/

// Examining convergence issues ....
//  with internal solver in IdealGases / Mixtures

  import Modelica.SIunits.*;
//  import ThermoS.Types.*;
  package Medium  = ThermoS.Media.MyGas;

//  constant Real AirComp[3] = {0.767,0.233,0.000};

  Medium.ThermodynamicState state(T(start=300),
                                  X(start=Medium.reference_X),
                                  p(start=1e5))   ;  // This is how you can have implicit start
  Temperature T ;
  SpecificEnthalpy h(start=1e5);
  parameter Integer s = 1 ; // works with s=0 , but fails with s=1 
// well ... with this version both work ...v1.11.0-dev.9+g4638cca
// import Medium = ThermoS.Media.MyGas is a nono ... you get weird behaviour. No simulation at all

  
equation
// Given Enthalpy find T of a Mix
     der(h) = 2.0e5 ;
     state.T = T ;

     if (s==1) then  // this does not converge (uses T_hX call
                     //   with internal SingleNonLinearEqn. solver
                     //  is this deficient???
        state = Medium.setState_phX(1e5, h, Medium.reference_X[1:2]);
     else          // but this does work .. no call to T_hX
       state.X = Medium.reference_X;
       state.p = 1e5;
       h = Medium.specificEnthalpy(state);
     end if;
end plant;
