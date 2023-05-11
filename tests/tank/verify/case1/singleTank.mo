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
  Feed    feed1(redeclare package Medium = Water); 
  Product    prod1(redeclare package Medium = Water); 
  Real U, T1 ;

equation

     connect (feed1.outlet, tank1.inlet) ;
     connect (tank1.outlet, prod1.inlet) ;


// Feed Stream Conditions
     feed1.mdot =  2.0/60;  
     feed1.T = 300 ; 
 
// Heat Addition and tank capacity
     tank1.Q_in = U * 0.5 * ((100+273) - tank1.Tf) ;
     tank1.tvol = 1000 ;

// Well we need to say where our inlet and outlets are  
     tank1.inlet.p = tank1.Pa ;  // This spec is need to use simple tank
     tank1.outlet.p = tank1.Pa ; // Ignoring static pressure effect (outlet can be anywhere
    
// Out flows of the tanks
     prod1.inlet.m_flow = feed1.mdot ;

// Reverse flow conditions of product streams
     prod1.T = 300 ; prod1.p = 2e5 ; 

     U = (4000 * 4.17) / 60 ;
     tank1.Pa = 1e5 ;  
     T1 = tank1.Tf - 273 ;
     

initial algorithm
    tank1.mf := 50 ;  tank1.Tf := 20 + 273 ; 
end plant;
