model plant
/*
  Author: Ravi Saripalli
  	 12th Feb 2015 
*/
  import ThermoS.Uops.*;
  import ThermoS.Types.*;
  import Water = Modelica.Media.Water.StandardWater;
  import ThermoS.Uops.Tanks.SimpleTank ;

  SimpleTank    tank1(redeclare package Medium = Water); 
  SimpleTank    tank2(redeclare package Medium = Water); 
  Feed    feed1(redeclare package Medium = Water); 
  Feed    feed2(redeclare package Medium = Water); 
  Product    prod1(redeclare package Medium = Water); 
  Product    prod2(redeclare package Medium = Water); 
  Real U, T1, T2 ;

equation

     connect (feed1.outlet, tank1.inlet) ;
     connect (tank1.outlet, prod1.inlet) ;

     connect (feed2.outlet, tank2.inlet) ;
     connect (tank2.outlet, prod2.inlet) ;

// Feed Stream Conditions
     feed1.mdot =  2.0/60;   feed2.mdot =  2.0/60;
     feed1.T = tank2.Tf ; feed1.Xi = tank2.outlet.Xi_outflow ;
     feed2.T = tank1.Tf ; feed2.Xi = tank1.outlet.Xi_outflow ;
 
// Heat Addition and tank capacity
     tank1.Q_in = U * 0.5 * ((100+273) - tank1.Tf) ;
     tank2.Q_in = U * 0.5 * ((120+273) - tank2.Tf) ;
     tank1.tvol = 1000 ; tank2.tvol = 1000 ;

// ===========================================================

/*These eqns should have been intrinsically added
  since mediaum state of SimpleTank has Pa pressure
  But I had to explicitly add these. What is even strange
  is that with 3.2.1 library only inlet port pressure linkage
  need to be specified
  But now with Lib 4.0.0 we need for all connected ports, inlet or outlet
*/

// Well we need to say where our inlet is 
     tank1.inlet.p = tank1.Pa ;  // This spec is need to use simple tank
     tank2.inlet.p = tank2.Pa ;  // This spec is need to use simple tank
    
// Don;t know why I got away without these two before
// Outlets need to know the pressure of the tank
     tank1.outlet.p = tank1.Pa ;
     tank2.outlet.p = tank2.Pa ;
// ===========================================================


// Out flows of the tanks
     prod1.inlet.m_flow = feed1.mdot ;
     prod2.inlet.m_flow = feed2.mdot ;

// Reverse flow conditions of product streams
     prod1.T = 300 ; prod1.p = 1e5 ; prod1.Xi = fill(0.0,0) ;
     prod2.T = 300 ; prod2.p = 1e5 ; prod2.Xi = fill(0.0,0) ;

     U = (4000 * 4.17) / 60 ;
     tank1.Pa = 1e5 ;  tank2.Pa = 1e5 ;
     T1 = tank1.Tf - 273 ;
     T2 = tank2.Tf - 273 ;


initial algorithm
    tank1.mf := 50 ;  tank1.Tf := 20 + 273 ; 
    tank2.mf := 100 ; tank2.Tf := 10 + 273 ; 
end plant;
