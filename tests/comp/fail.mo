model plant
  /*
  Author: Ravi Saripalli
  	11st Oct 2025 
    Compressor Tests 
 atm -> comp -> tank -> vout -> atm
But this model fails .. beats me. 
It runs if I use Media with reducedX = false
Definite bug in the 
*/
  import ThermoS.Types.* ;
  import ThermoS.Media.MyGas ;
  import ThermoS.Uops.Valves.* ;
  import ThermoS.Uops.Machines.FlowMachine ;
  import ThermoS.Uops.Tanks.*;
  import Modelica.Units.SI.MassFlowRate;

  
  constant    Real CV = 200 * (1e-3/60) / sqrt (10000) ;
  constant    Real Air[MyGas.nXi] = {0.7, 0.2, 0.1} ; // Changed to full Xi
  Reservoir   atm (redeclare package Medium = MyGas,
                         p = 1e5, T = 300, Xi = Air);
  Valve       vout (redeclare package Medium = MyGas, cv = CV);
  GasTank     tank (redeclare package Medium = MyGas, vol = 0.2); 
  FlowMachine comp (redeclare package Medium = MyGas, isComp = true, eff = 1) ;
  constant    MassFlowRate  maxFlow = 200 * 1e-3 / 60 ;

equation
     connect (atm.port, comp.inlet) ;
     connect (comp.outlet, tank.inlet);
     connect (tank.outlet, vout.inlet) ;
     connect (vout.outlet, atm.port) ;
     vout.po = 100 ;
     // Compressor Curve (just linear for now
     comp.outlet.p = ( 1.9 - (comp.inlet.m_flow / maxFlow) ) * comp.inlet.p ; 
     // Assuming 2m2 surface area, 15 outside heattransfer coeff
     tank.Q_loss = 15 * 2 * (tank.T - (15 + 273)) ;

initial algorithm
    tank.T := 300 ;  // Initial Temperature
    tank.p := 1.0e5  ;  // Initial Pressure
    tank.Xi := Air ;
end plant;
