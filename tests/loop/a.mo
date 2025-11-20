model plant
/*
  Author: Ravi Saripalli
  	 5th Nov. 2025 
   Test Rig Pressurization loop 

   Big lesson:  An unconnected fluid port would always
                have zero flow through it. So one cant
                arbitrarily set its m_flow. I thought
                an unconnected port is unspecified port.
                But as per Modelica it is a closed port.

uptoCMp  + TeeCmp to atm
*/
  import ThermoS.Types.* ;
  import Modelica.Units.SI.* ;
  import ThermoS.Media.* ;
  import ThermoS.Uops.Valves.* ;
  import ThermoS.Uops.Tanks.* ;
  import ThermoS.Uops.Machines.* ;
  import ThermoS.Uops.HeatExch.HeatX ;

 // import ThermoS.Uops.Controller ;

  constant    Real Air[MyAir.nXi] = {0.7} ;
  // A 25mm valve at 200 lpm gives about 100Pa pressure drop
  //   when half open.
  constant    Real CV  = (200 * 1e-3/60 ) / sqrt (100) ;
  constant    MassFlowRate  maxFlow = 1 ;
  constant    Real maxPr = 1.9 ;

  // Ambient Conditions
  Reservoir   atm (redeclare package Medium = MyAir,
                         p = 1e5, T = 300, Xi = Air);
  GasTank     volUp (redeclare package Medium = MyAir, vol = 0.2 ); 
  RealValve   vUp   (redeclare package Medium = MyAir , cv = CV);
  RealValve   vDown   (redeclare package Medium = MyAir , cv = CV);
  RealValve   vBlow (redeclare package Medium = MyAir , cv = CV);
  FlowMachine  cmpTrim (redeclare package Medium = MyAir, 
                    isComp = true, eff = 0.95) ;  
  HeatX  hx (redeclare package Medium = MyAir, cf = 1.0e-3, 
                       A_wf = 1,  h_wf = 150, 
                       w_m = 1, w_cp = 420, holdup = 100e-3);
  portMixer TeeBlow (redeclare package Medium = MyAir, 
                        vol=1e-3, N=3, Adiabatic = true) ;
  portMixer TeeHx (redeclare package Medium = MyAir, 
                        vol=1e-3, N=3, Adiabatic = true) ;

  portMixer TeeComp (redeclare package Medium = MyAir, 
                         N=2, Adiabatic = true) ;
equation
   
     // Trimming Compressor + blow off valve
     connect (atm.port, cmpTrim.inlet) ;
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

     vBlow.spo = 90 ;
     vUp.spo = 90 ;
     vDown.spo = 90 ;

     cmpTrim.Ws = max (0, 100 * (1.7 - volUp.p/1e5)/1.7) ;
     cmpTrim.outlet.p =  (maxPr  - (maxPr - 1)*(cmpTrim.inlet.m_flow/maxFlow)) * cmpTrim.inlet.p ; 
     hx.Q_ew = -10 ; // External to wall heat transfer
     volUp.Q_loss = 30 ;

initial algorithm
    volUp.T := 300 ;  // Initial Temperature
    volUp.p := 1e5  ;  // Initial Pressure
    volUp.Xi := Air ;

    hx.Tf := 300 ;  // Heater fluid Temp
    hx.Tw := 300 ;  // Heater wall Temp
    
     vBlow.po := 5 ;
     vUp.po := 5 ;
     vDown.po := 5 ;

     TeeBlow.p := 1e5 ;
     TeeBlow.T := 300 ;
     TeeBlow.Xi := Air ;

     TeeHx.p := 1e5 ;
     TeeHx.T := 300 ;
     TeeHx.Xi := Air ;

     TeeComp.p := 1e5 ;
     TeeComp.T := 300 ;
     TeeComp.Xi := Air ;
end plant;
