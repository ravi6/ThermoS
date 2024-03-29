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

// Set pressure in the tank at ports
   tank1.inlet.p = tank1.Pa ;
   tank2.inlet.p = tank1.Pa ;
   tank1.outlet.p = tank1.Pa ;
   tank2.outlet.p = tank1.Pa ;


initial algorithm
    tank1.mf := 50 ;  tank1.Tf := 20 + 273 ; 
    tank2.mf := 100 ; tank2.Tf := 10 + 273 ; 
end plant;
