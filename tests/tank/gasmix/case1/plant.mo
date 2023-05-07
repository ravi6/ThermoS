model plant
/*
  Author: Ravi Saripalli
  	 12th Feb 2015 
*/
  import ThermoS.Types.*;
  import ThermoS.Uops.*;
  import ThermoS.Media.MyGas;
  import ThermoS.Uops.Valves.Valve;
  import ThermoS.Uops.Tanks.OpenTank ;

  constant    Real Air[MyGas.nXi] = {0.79, 0.21} ;
  OpenTank    tank(redeclare package Medium = MyGas, in_pos = 0.44); 
  Reservoir   lake(redeclare package Medium = MyGas, p=1e5, T=300, Xi=Air); // Reservoir 1
  Feed        supply (redeclare package Medium = MyGas); // InletFlow to tank
  Valve       valve (redeclare package Medium = MyGas) ; //, cv = 1e-3 / sqrt(100)) ;

equation
     supply.mdot = 5 + 4.95 * sin(6*time) ; 
     supply.T = 300 ; // for a force feed you need flow, temp and comp
     supply.Xi = fill(1.0/MyGas.nS, MyGas.nXi) ;

     connect (supply.outlet, tank.inlet) ;
     connect (tank.outlet, valve.inlet) ;
     connect (valve.outlet, lake.port) ;

     tank.hcoef = 150 ;
     tank.Pa   = 1e5 ;  tank.Ta =300 ;
     valve.po  = 50 ;

initial algorithm
    tank.pFull := 50 ; tank.Tf := 400 ; 

end plant;
