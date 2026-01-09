model plant
/*
  Author: Ravi Saripalli
  	 5th Nov. 2025 
   Test Rig Pressurization loop 
   Upto Heat Exchager outlet
*/

 // import ThermoS.Uops.Controller ;
  import ThermoS.Types.* ;
  import Modelica.Units.SI.* ;
  import ThermoS.Media.* ;
  import ThermoS.Uops.Valves.* ;
  import ThermoS.Uops.Tanks.* ;
  import ThermoS.Uops.Machines.* ;
  import ThermoS.Uops.HeatExch.HeatX ;

  constant    Real Air[MyAir.nXi] = {0.7} ;
  constant    Real CV  = (200 * 1e-3/60 ) / sqrt (1000) ;
  constant    MassFlowRate  maxFlow = 1 ;

  // Ambient Conditions
  Reservoir   atm (redeclare package Medium = MyAir,
                         p = 1e5, T = 300, Xi = Air);

  // Depressurizing blow off valve ? 
  RealValve   vBlow (redeclare package Medium = MyAir , cv = CV);

  // Pressure controlling Compressor
  FlowMachine  cmpTrim (redeclare package Medium = MyAir, 
                    isComp = true, eff = 0.95) ;  

  // Cooler after compression
  HeatX hx (redeclare package Medium = MyAir, cf = 1.0e-3, 
                       A_wf = 1,  h_wf = 150, 
                       w_m = 1, w_cp = 420, holdup = 10e-3);

  // Blow valve Tee junction
  portMixer TeeBlow (redeclare package Medium = MyAir, 
                        vol=1000e-3, N=3, Adiabatic = true) ;

equation
   
     // Trimming Compressor + blow off valve
     connect (atm.port, cmpTrim.inlet) ;
     connect (cmpTrim.outlet, TeeBlow.port[1]) ;
     connect (TeeBlow.port[2], vBlow.inlet) ;
     connect (vBlow.outlet, atm.port) ; 
 
     // Cooler before entering test section
     connect (TeeBlow.port[3], hx.inlet) ; 
     connect (hx.outlet, atm.port)  ;

     vBlow.spo = 90 ;
     cmpTrim.Ws = 100;
     cmpTrim.outlet.p =  (1.2  - (1.2 - 1) *(cmpTrim.inlet.m_flow/maxFlow)) * cmpTrim.inlet.p ; 

     hx.Q_ew = -1000 ; // External to wall heat transfer

initial algorithm
    hx.Tf := 300 ;  // Heater fluid Temp
    hx.Tw := 300 ;  // Heater wall Temp
    vBlow.po := 10 ;
    TeeBlow.p := 1e5 ;
    TeeBlow.T := 300 ;
    TeeBlow.medium.Xi := Air ;

/*
    volUp.T := 300 ;  // Initial Temperature
    volUp.p := 0.999e5  ;  // Initial Pressure
    volUp.Xi := Air ;
    volDown.T := 300 ;  // Initial Temperature
    volDown.p := 0.999e5  ;  // Initial Pressure
    volDown.Xi := Air ;
    
     vBlow.po := 5 ;
     vUp.po := 5 ;
     vDown.po := 5 ;
     vBack.po := 5 ;

 */   
end plant;
