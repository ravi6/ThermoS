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
  Reservoir   suction  (redeclare package Medium = MyGas, p=1e5, T=300, Xi=Air);
  Reservoir   atm (redeclare package Medium = MyGas, p=1e5, T=300, Xi=Air);
  Valve       valve (redeclare package Medium = MyGas , 
                     cv = (1500e-3/60) / sqrt(4e5)) ;
  CompressorBasic  comp (redeclare package Medium = MyGas, 
                                 mdot  =  1.0 , pr = 2,  n = 1.4) ;
  GasTank     tank (redeclare package Medium = MyGas, vol = 0.2 , Q_in=0); 

equation
     connect (suction.port, comp.inlet) ;
     connect (comp.outlet, tank.inlet) ;
     connect (tank.outlet, valve.inlet) ;
     connect (valve.outlet, atm.port) ;
     
initial algorithm
    tank.T := 300 ;  // Initial Temperature
    tank.p := 1e5 ;  // Initial Pressure
    tank.Xi := Air ;
    comp.inlet.m_flow := 0 ;
end plant;
