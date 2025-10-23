model plant
/*
  Author: Ravi Saripalli
  	11st Oct 2025 
    Compressor Tests 
*/
  import ThermoS.Types.* ;
  import ThermoS.Media.MyGas ;
  import ThermoS.Uops.Tanks.GasTank ;
  import ThermoS.Uops.Valves.Valve ;
  import ThermoS.Uops.CompressorBasic ;
  import ThermoS.Uops.Reservoir ;

  constant    Real Air[MyGas.nXi] = {0.7, 0.2} ;
  Reservoir   suction  (redeclare package Medium = MyGas,
                         p = 1e5, T = 300, Xi = Air);
  Reservoir   atm (redeclare package Medium = MyGas,
                         p = 1e5, T = 300, Xi = Air);
  Valve       vin (redeclare package Medium = MyGas , 
                     cv = (1500e-3/60) / sqrt(4e5)) ;
  Valve       vout (redeclare package Medium = MyGas , 
                     cv = (1500e-3/60) / sqrt(4e5)) ;
  GasTank     tank (redeclare package Medium = MyGas, 
                            vol = 0.2 , Q_in=0); 
  CompressorBasic  comp (redeclare package Medium = MyGas, 
                              holdup = 0.001, eff = 0.95) ;

equation
     connect (suction.port, comp.inlet) ;
     connect (comp.outlet, vin.inlet) ;
     connect (vin.outlet, tank.inlet) ;
     connect (tank.outlet, vout.inlet) ;
     connect (vout.outlet, atm.port) ;
     vin.po = 50 ;
     vout.po = 50 ;
     comp.Ws = -100;
    comp.inlet.m_flow = 1 ;

initial equation
    tank.T = 300 ;  // Initial Temperature
    tank.p = 1.0e5  ;  // Initial Pressure
    tank.Xi = Air ;
    comp.state_out.T = 300 ;
    comp.state_out.p = 1.0e5 ;
//    comp.U = comp.Medium.specificInternalEnergy(comp.state_out) ;
end plant;
