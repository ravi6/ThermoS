model plant
/*
  Author: Ravi Saripalli
  	 12th Feb 2015
*/
  import ThermoS.Types.*;
  import ThermoS.Uops.*;
  import ThermoS.Media.JP8;
  import ThermoS.Uops.Tanks.OpenTank ;
  import ThermoS.Uops.Valves.Valve ;
  import ThermoS.Uops.Reservoir;
  import ThermoS.Uops.PumpBasic ;

  OpenTank    tank(redeclare package Medium = JP8, in_pos = 0.44); 
  Reservoir   lake(redeclare package Medium = JP8, p= 1e5, T = 300); // Xi will be null for single component fluid
                                                    // but you need to specify for multicomponents Xi=fill(1./JP8.nS, JP8.nXi)); 
  Reservoir   dam(redeclare package Medium = JP8, p = 2e5, T = 300); 
  Valve       vout (redeclare package Medium = JP8) ; //, cv = 1e-3 / sqrt(100)) ;
  Valve       vin (redeclare package Medium = JP8, cv = 0.5e-2 / sqrt(100)) ;
  PumpBasic   pump (redeclare package Medium = JP8); 

  Real dpvh ;

equation
     connect (dam.port, vin.inlet) ;
     connect (vin.outlet, tank.inlet) ;
     connect (tank.outlet, vout.inlet) ;
     connect (vout.outlet, pump.inlet) ;
     connect (pump.outlet, lake.port) ;

     tank.hcoef = 150 ;
     tank.Pa   = 1e5 ;  tank.Ta =300 ;
     vin.po = 80 ; vout.po = 50 ;
     pump.Q = 0 ; pump.Ws = 20 ;
     dpvh = -(vout.outlet.p - vout.inlet.p) / (pump.Medium.density(pump.state) * 9.81) ;
initial algorithm
    tank.pFull := 80 ; tank.Tf    := 400 ; 
    pump.T := 300 ;
end plant;
