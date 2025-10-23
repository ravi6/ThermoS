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

  constant    Real Air[MyGas.nXi] = {0.7, 0.2} ;
  Reservoir   suction  (redeclare package Medium = MyGas,
                               p = 1e5, T = 300, Xi = Air);
  Reservoir   atm (redeclare package Medium = MyGas,
                               p = 3e5, T = 300, Xi = Air);
  CompressorBasic  comp (redeclare package Medium = MyGas, 
                              holdup = 0.001, eff = 0.95) ;

equation
     connect (suction.port, comp.inlet) ;
     connect (comp.outlet, atm.port) ;
     comp.Ws = -100;

initial equation
//    comp.U = comp.Medium.specificInternalEnergy(comp.state_out) ;
end plant;
