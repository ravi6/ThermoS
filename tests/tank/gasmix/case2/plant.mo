model plant
/*
  Author: Ravi Saripalli
  	 12th Feb 2015
*/
  import ThermoS.Types.*;
  import ThermoS.Uops.*;
  import ThermoS.Media.MyGas;
  import ThermoS.Uops.Tanks.OpenTank ;
  import ThermoS.Uops.Valves.Valve ;
  import ThermoS.Uops.Reservoir;

  constant Real Air[MyGas.nXi]={0.79, 0.21} ; 
  OpenTank    tank(redeclare package Medium = MyGas, in_pos = 0.44); 
  Reservoir   lake(redeclare package Medium = MyGas, p = 1e5, T = 300, Xi=Air); // Xi will be null for single component fluid
                                                    // but you need to specify for multicomponents Xi=fill(1./MyGas.nS, MyGas.nXi)); 
  Reservoir   dam(redeclare package Medium = MyGas, p = 2e5, T = 300, Xi=Air); 
  Valve       vout (redeclare package Medium = MyGas) ; //, cv = 1e-3 / sqrt(100)) ;
  Valve       vin (redeclare package Medium = MyGas, cv = 0.5e-2 / sqrt(100)) ;

equation
     connect (dam.port, vin.inlet) ;
     connect (vin.outlet, tank.inlet) ;
     connect (tank.outlet, vout.inlet) ;
     connect (vout.outlet, lake.port) ;
     tank.hcoef = 150 ;
     tank.Pa   = 1e5 ;  tank.Ta =300 ;
     vin.po = 100 ; vout.po = 100 ;

initial algorithm
    tank.pFull := 80 ; tank.Tf    := 400 ; 
end plant;
