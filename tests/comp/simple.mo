model plant
/*
  Author: Ravi Saripalli
  	11st Oct 2025 
    Compressor Tests 
*/
  import ThermoS.Types.* ;
  import ThermoS.Media.MyGas ;
  import ThermoS.Uops.CompressorBasic ;
  import ThermoS.Uops.Reservoir ;

  constant    MyGas.MassFraction  Air[MyGas.nXi] = {0.7, 0.2} ;
  Reservoir   suction  (redeclare package Medium = MyGas,
                         p = 1e5, T = 300, Xi = Air);
  Reservoir   disch  (redeclare package Medium = MyGas,
                         p = 2e5, T = 300, Xi = Air);
  CompressorBasic  comp (redeclare package Medium = MyGas, 
                             eff = 0.95) ;

equation
     connect (suction.port, comp.inlet) ;
     connect (comp.outlet, disch.port) ;
     comp.Ws = 1000 ;

initial equation
end plant;
