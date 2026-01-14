model plant
/*
  Author: Ravi Saripalli
  	11st Oct 2025 
    Compressor Tests 
 atm -> comp -> htr -> tank -> vout -> atm
*/
  import ThermoS.Types.* ;
  import ThermoS.Media.MyGas ;
  import ThermoS.Uops.Valves.* ;
  import ThermoS.Uops.Machines.FlowMachine ;
  import ThermoS.Uops.Tanks.*;
  import ThermoS.Uops.HeatExch.HeatX;
  import Modelica.Units.SI.MassFlowRate;

  
  constant    MassFlowRate  maxFlow = 200 * 1e-3 / 60 ;
  constant    Real CV = 200 * (1e-3/60) / sqrt (10000) ;
  constant    Real Air[2] = {0.7, 0.2} ;

  Reservoir   atm (redeclare package Medium = MyGas,
                         p = 1e5, T = 300, Xi = Air);
  Valve       vout (redeclare package Medium = MyGas, cv = CV);
  GasTank     tank (redeclare package Medium = MyGas, vol = 0.2); 
  FlowMachine comp (redeclare package Medium = MyGas, isComp = true, eff = 1) ;
  HeatX htr (redeclare package Medium = MyGas, cf = 1.0e-3, 
                                A_wf = 1,  h_wf = 150, 
                                w_m = 1, w_cp = 420, holdup = 50);

equation
     connect (atm.port, comp.inlet) ;
     connect (comp.outlet, htr.inlet);
     connect (htr.outlet, tank.inlet) ;
     connect (tank.outlet, vout.inlet) ;
     connect (vout.outlet, atm.port) ;
     vout.po = 100 ;
     comp.Ws = 1000 ;
     // Assuming 2m2 surface area, 15 outside heattransfer coeff
     tank.Q_loss = 15 * 2 * (tank.T - (15 + 273)) ;
     htr.Q_ew = 500 ;

initial algorithm
    tank.T := 300 ;  // Initial Temperature
    tank.p := 1.0e5  ;  // Initial Pressure
    tank.Xi := {0.7, 0.2} ; 
    htr.Tf :=300 ;
    htr.Tw :=300 ; 
    htr.medium.p := 1e5;
end plant;
