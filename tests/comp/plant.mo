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
  import ThermoS.Uops.FlowMachine ;
  import ThermoS.Uops.Reservoir ;

  constant    Real Air[MyGas.nXi] = {0.7, 0.2} ;
  Reservoir   atm (redeclare package Medium = MyGas,
                         p = 1e5, T = 300, Xi = Air);
  Valve       vin (redeclare package Medium = MyGas, 
                     cv = (1000/60) / sqrt(.01e5)) ;
  Valve       vout (redeclare package Medium = MyGas, 
                     cv = (1000/60) / sqrt(.01e5)) ;
  GasTank     tank (redeclare package Medium = MyGas, 
                            vol = 0.2 , Q_in=0); 
  FlowMachine comp (redeclare package Medium = MyGas, 
                       isComp = true, prat = 1.2, eff = 1) ;

equation
     connect (atm.port, comp.inlet) ;
     connect (comp.outlet, vin.inlet) ;
     connect (vin.outlet, tank.inlet) ;
     connect (tank.outlet, vout.inlet) ;
     connect (vout.outlet, atm.port) ;
     vin.po = 50 ;
     vout.po = 50 ;
     comp.Ws = 100 ;

initial algorithm
    tank.T := 300 ;  // Initial Temperature
    tank.p := 1.0e5  ;  // Initial Pressure
    tank.Xi := Air ;
end plant;
