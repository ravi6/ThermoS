model plant
/*
  Author: Ravi Saripalli
  	11st Oct 2025 
    Compressor Tests 
 atm -> comp -> vin - tank - vout - atm
*/
  import ThermoS.Types.* ;
  import ThermoS.Media.MyGas ;
  import ThermoS.Uops.Valves.Valve ;
  import ThermoS.Uops.Machines.FlowMachine ;
  import ThermoS.Uops.Tanks.*;
  import Modelica.Units.SI.MassFlowRate;

  
  constant    Real CV = 200 * (1e-3/60) / sqrt (10000) ;
  constant    Real Air[MyGas.nXi] = {0.7, 0.2} ;
  Reservoir   atm (redeclare package Medium = MyGas,
                         p = 1e5, T = 300, Xi = Air);
  Valve       vin  (redeclare package Medium = MyGas, cv = CV);
  Valve       vout (redeclare package Medium = MyGas, cv = CV);
  GasTank     tank (redeclare package Medium = MyGas, vol = 0.2); 
  FlowMachine comp (redeclare package Medium = MyGas, isComp = true, eff = 1) ;
  constant    MassFlowRate  maxFlow = 200 * 1e-3 / 60 ;
  constant    Real          prMax = 1.01 ;
  MassFlowRate mflow (start=0.0001) ;

equation
     connect (atm.port, comp.inlet) ;
     connect (comp.outlet, vin.inlet) ;
     connect (vin.outlet, tank.inlet) ;
     connect (tank.outlet, vout.inlet) ;
     connect (vout.outlet, atm.port) ;
     vin.po = 100 ;
     vout.po = 100 ;

     // Compressor Curve (just linear for now
     der(mflow - comp.inlet.m_flow) = - mflow + comp.inlet.m_flow ;
     comp.outlet.p = ( prMax - (prMax - 1) * (mflow / maxFlow) ) * comp.inlet.p   ; 
     comp.Ws = 200 * (1 - exp(- time)) ;
     // Assuming 2m2 surface area, 15 outside heattransfer coeff
     tank.Q_loss = 15 * 2 * (tank.T - (15 + 273)) ;

initial equation
    tank.medium.state = MyGas.setState_pTX (1e5, 300, Air) ;
    comp.medium.state = MyGas.setState_pTX (1e5, 300, Air) ;
    comp.inlet.m_flow = 0.001 ;
    // tank.T = 300 ;  // Initial Temperature
    // tank.p = 1.0e5  ;  // Initial Pressure
    // tank.Xi = Air ;
    // comp.medium.T = 300 ;
    // comp.medium.p =  1e5 ;
    // comp.medium.Xi = Air ;
end plant;
