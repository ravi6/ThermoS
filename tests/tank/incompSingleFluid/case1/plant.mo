model plant
/*
  Author: Ravi Saripalli
  	 12th Feb 2015 
*/
  import ThermoS.Uops.*;
  import ThermoS.Uops.Valves.*;
  import ThermoS.Types.*;
  import ThermoS.Media.JP8;
  import ThermoS.Uops.Tanks.OpenTank ;

  OpenTank    tank(redeclare package Medium = JP8, in_pos = 0.44); 
  Reservoir   lake(redeclare package Medium = JP8, p=1e5, T=300); // Reservoir 1
  Feed        supply (redeclare package Medium = JP8); // InletFlow to tank
  Valve       valve (redeclare package Medium = JP8) ; //, cv = 1e-3 / sqrt(100)) ;

equation
     supply.mdot = 3 + 0.2 * sin(0.1*time) ; supply.T = 300 ; // for a force feed you need flow, temp and comp
                                                            // composition in this case non-existent
                                                           // nS=1 and reducedX=true

     connect (supply.outlet, tank.inlet) ;
     connect (tank.outlet, valve.inlet) ;
     connect (valve.outlet, lake.port) ;

     tank.hcoef = 150 ;
     tank.Pa    = 1e5 ;  tank.Ta = 300 ;
     valve.po   = 60 ;

initial algorithm
    tank.pFull := 50 ; tank.Tf := 400 ; 

end plant;
