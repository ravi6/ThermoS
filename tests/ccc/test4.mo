model plant
/*
  Author: Ravi Saripalli
  	 5th Nov. 2025 
   Test Rig Pressurization loop 
   
   Test compressor omitted.


   Big lesson:  An unconnected fluid port would always
                have zero flow through it. So one cant
                arbitrarily set its m_flow. I thought
                an unconnected port is unspecified port.
                But as per Modelica it is a closed port.

uptoCMp  + TeeCmp to atm
*/
  import Rig.* ;

 // import ThermoS.Uops.Controller ;

  // A 25mm valve at 200 lpm gives about 100Pa pressure drop
  //   when half open.
  constant    Real CV  = (200 * 1e-3/60 ) / sqrt (100) ;
  constant    Real  maxFlow = 1 ;
  constant    Real maxPr = 1.9 ;

  // Ambient Conditions
  Res     atm (redeclare package Medium = Air, p = 1e5, T = 300);
  Valve   vin    (redeclare package Medium = Air , cv = CV);
  Tank    volUp (redeclare package Medium = Air, vol = 0.2 ); 
  Valve   vUp   (redeclare package Medium = Air , cv = CV);
  Valve   vDown   (redeclare package Medium = Air , cv = CV);
  Valve   vBlow (redeclare package Medium = Air , cv = CV);
  Comp    cmpTrim (redeclare package Medium = Air, eff = 0.95) ;  
  HeatX  hx (redeclare package Medium = Air, cf = 1.0e-3, 
                       A_wf = 1,  h_wf = 150, 
                       w_m = 1, w_cp = 420, holdup = 100e-3);
  Tee TeeBlow (redeclare package Medium = Air) ;
  Tee TeeHx (redeclare package Medium = Air) ;
  Tee TeeComp (redeclare package Medium = Air, N=2) ;

equation
   
     // Trimming Compressor + blow off valve
     connect (atm.port, vin.inlet);
     connect (vin.outlet, cmpTrim.inlet) ;
     connect (cmpTrim.outlet, TeeBlow.port[1]) ;
     connect (TeeBlow.port[2], vBlow.inlet) ;
     connect (vBlow.outlet, atm.port) ; 
 
     // Cooler before entering test section
     connect (TeeBlow.port[3], hx.inlet) ; 
     connect (hx.outlet, TeeHx.port[1])  ;

     // Flow split between Up and Downstream test section
     connect (TeeHx.port[2], vUp.inlet) ;
     connect (TeeHx.port[3], vDown.inlet) ;
    
     // Upstream Test chamber connections
     connect (vUp.outlet, volUp.inlet);
     connect (volUp.outlet, atm.port);
     connect (vDown.outlet, TeeComp.port[1]) ;
     connect (TeeComp.port[2], atm.port) ;
//     connect (TeeComp.port[3], atm.port) ;

     vin.po = 90 ;
     vBlow.po = 10 ;
     vUp.po = 90 ;
     vDown.po = 90 ;

    cmpTrim.pr = 1.1 ;
//     cmpTrim.outlet.p =  (maxPr  - (maxPr - 1)*(cmpTrim.inlet.m_flow/maxFlow)) * cmpTrim.inlet.p ; 

     hx.Q_ew = 10 ; // External to wall heat transfer
     volUp.Q_loss = 30 ;

initial algorithm
    volUp.T := 300 ;  // Initial Temperature
    volUp.p := 1e5  ;  // Initial Pressure

    hx.Tf := 300 ;  // Heater fluid Temp
    hx.Tw := 300 ;  // Heater wall Temp
    hx.p  := 1e5 ; 

end plant;
