model plant
/*
  Author: Ravi Saripalli
  	 12th Feb 2015
*/
  import ThermoS.Types.*;
  import ThermoS.Uops.*;
  import ThermoS.Media.JP8;
  import ThermoS.Uops.Tanks.OpenTank ;
  import Valve = ThermoS.Uops.Valves.NonReturnValve ;
  import ThermoS.Uops.Reservoir;

  OpenTank     tank(redeclare package Medium = JP8, in_pos = 0.3, h=1.0, d=0.1 ); 
  Valve        valve (redeclare package Medium = JP8, cv = 1e-4 / sqrt(100)) ;
  PumpBasic    pump (redeclare package Medium = JP8, holdup = 1e-4); 
  HeaterCooler htr(redeclare package Medium = JP8, cf = 1e-3, 
                                        A_wf = 1,  h_wf = 150, 
                                        w_m = 1, w_cp = 420, holdup = 1e-4);
  Reservoir    lake(redeclare package Medium = JP8, p = 1e5, T = 300) ;

  record recdelp
      Real pump ;
      Real valve ;
      Real heater ;
  end recdelp;

  recdelp delp ;


/*
                --------<-------------+                  
               |                      |
             | v    |                |>| 
             | !    |                |>| 
             \ Tank /                |>| Hx
              \ _ _/                 |>|
                 |                    |    @
                 !---->--(*)----------+---)!(----> ~~`~~
                         /_\                        ~~  Lake
                        pump             valve
*/

equation
     connect (tank.outlet, pump.inlet) ;
     connect (pump.outlet, htr.inlet) ;
     connect (htr.outlet, tank.inlet) ;
     connect (pump.outlet, valve.inlet) ;
     connect (valve.outlet, lake.port) ;
     

     tank.hcoef = 150 ; tank.Pa   = 1e5 ;  tank.Ta =300 ;
 //    valve.po = 80 ;
    pump.Q = 0 ; pump.head = 1.2  * (1-exp(-time/120)) ;
     htr.Q_ew =  2000 * (1 - exp(-time/40)) ;

     delp.pump  = (pump.outlet.p - pump.inlet.p) / 1e3 ;
     delp.valve  = (valve.inlet.p - valve.outlet.p) / 1e3 ;
     delp.heater  = (htr.inlet.p - htr.outlet.p) / 1e3 ;

initial algorithm
    tank.pFull := 80 ; tank.Tf    := 400 ; 
    htr.Tf :=300 ; htr.Tw :=300 ; 
    pump.T := 300 ;  
end plant;
