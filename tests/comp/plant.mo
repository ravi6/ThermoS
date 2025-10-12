model plant
  import ThermoS.Media.MyGas;
  import ThermoS.Uops.CompressorBasic;

  constant Real Air[MyGas.nXi] = {0.7, 0.2};

  CompressorBasic comp(
    redeclare package Medium = MyGas,
    pr = 2,
    n  = 1.4
  );

  // Dynamic inlet pressure
  Real p_target = 1e5;       // target inlet pressure [Pa]
  Real tau = 0.1;             // time constant for inlet pressure dynamics
  Real x;                     // example independent dynamic state

equation
  // First-order dynamic for inlet pressure
  der(comp.inlet.p) = (p_target - comp.inlet.p)/tau;

  // Feed the compressor with enthalpy and composition based on current inlet pressure
  comp.inlet.Xi_outflow = Air;
  comp.inlet.h_outflow = MyGas.specificEnthalpy(
    MyGas.setState_pTX(comp.inlet.p, 300, Air)
  );

  // Example dummy ODE
  der(x) = 1;

  // Output to prevent pruning
  x = x; // trivial, ensures at least one output

initial algorithm
  comp.inlet.p := 1e5; // start near target
  x := 0;

end plant;
