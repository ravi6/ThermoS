model plant
/*
  Author: Ravi Saripalli
Get Temperature of gas mixture given enthalpy, pressure and composition

*/
    
    
  // package Medium = Modelica.Media.IdealGases.MixtureGases.CombustionAir(X(start={0.79,0.21}), T(start=300), p(start=1e5));
  package Medium = ThermoS.Media.MyAir ;
  import Modelica.Media.IdealGases.Common.MixtureGasNasa.h_TX;

  Medium.ThermodynamicState state(T(start=300), X(start={0.79, 0.21}), p(start=1e5)) ;
  Medium.SpecificEnthalpy h ;
    
  constant Integer option = 3  ;

equation

    h = 600000;

if option == 1 then     // works 
// invariably calls T_hX function which in turn uses internal singleNonlinear eqn. solver
//  backend Translation of which fails on some  instances
    state = Medium.setState_phX(1e5, h, {0.79, 0.21});
elseif option == 2 then  // works ... and avoids T_hX calls
    state.p = 1e5 ;
    state.X = {0.79, 0.21} ;
    h  -  Medium.specificEnthalpy(state) = 0;
else 
     state.p = 1e5 ;
     state.X = {0.79, 0.21} ;
      h  -  h_TX(state.T, state.X) = 0 ; // fails why
    //  h * state.p * state.T - state.T * state.p * h_TX(state.T, state.X) = 0 ;  // results in crap result of state.T=0
    //  h *  state.T - state.T *  h_TX(state.T, state.X) = 0 ;  // results in segmentation fault
 end if;                                     
end plant;
