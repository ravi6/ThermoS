model plant
/*
  Author: Ravi Saripalli
  	11st Oct 2025 
    Compressor Tests 
*/
  import ThermoS.Types.* ;
  import ThermoS.Media.MyGas ;
  import ThermoS.Uops.Feed ;
//  import ThermoS.Uops.Product ;
  import ThermoS.Uops.Tanks.GasTank ;
  import ThermoS.Uops.Valves.Valve ;
  //import ThermoS.Uops.CompressorBasic ;
  import ThermoS.Uops.Reservoir ;

  constant    Real Air[MyGas.nXi] = {0.7, 0.2} ;
  Feed      suction (redeclare package Medium = MyGas);
  Valve     valve (redeclare package Medium = MyGas , 
                     cv = (1500e-3/60) / sqrt(4e5)) ;
 // CompressorBasic  comp (redeclare package Medium = MyGas, 
                         //         pr = 2,  n=1.4) ;
  GasTank     tank (redeclare package Medium = MyGas, vol = 0.2 , Q_in=0); 
  Reservoir   atm (redeclare package Medium = MyGas, p=1e5, T=300, Xi=Air);
  constant Real alpha = 1.0e-2 ;
  Real x (start=0, min = 0 , max = 1);

equation
     connect (suction.outlet, tank.inlet) ;
//     connect (tank.outlet, comp.inlet) ;
//     connect (comp.outlet, valve.inlet) ;
     
     connect (tank.outlet, valve.inlet);
     connect (valve.outlet, atm.port) ;
     suction.T = 300  ;
     suction.Xi = Air ;
     valve.po = 50 ;
     tank.inlet.m_flow =  x ;
     x = alpha * tank.p ;

initial algorithm
    tank.T := 300 ;  // Initial Temperature
    tank.p := 1e5 ;  // Initial Pressure
    tank.Xi := Air ;
end plant;
