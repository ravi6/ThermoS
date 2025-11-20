model plant
/*
  Author: Ravi Saripalli
  	 5th Nov. 2025 
   Test Rig Pressurization loop 
*/

  import ThermoS.Types.* ;
  import Modelica.Units.SI.* ;
  import ThermoS.Media.* ;
  import ThermoS.Uops.Valves.* ;
  import ThermoS.Uops.Tanks.* ;
  import ThermoS.Uops.Machines.* ;
  import ThermoS.Uops.HeatExch.HeatX ;

  constant    Real Air[MyAir.nXi] = {0.7} ;
  constant    Real CV  = (200 * 1e-3/60 ) / sqrt (100) ;
  constant    MassFlowRate  maxFlow = 1 ;
  constant    Real maxPr = 1.9 ;

  // Ambient Conditions
  Reservoir   atm (redeclare package Medium = MyAir,
                         p = 1e5, T = 300, Xi = Air);

  // Equivalent Test champer Volumes (across test Compressor)
  GasTank     volUp (redeclare package Medium = MyAir, vol = 0.2 ); 
  GasTank     volDown (redeclare package Medium = MyAir, vol = 0.2 ); 

  // Control valves feeding upstream and downstream 
  // holdup volumes across test compressor
  RealValve   vUp   (redeclare package Medium = MyAir , cv = CV);
  RealValve   vDown (redeclare package Medium = MyAir , cv = CV);

  // Back pressure control valve
  RealValve   vBack (redeclare package Medium = MyAir , cv = CV);

  // Depressurizing blow off valve ? 
  RealValve   vBlow (redeclare package Medium = MyAir , cv = CV);

  // Compressor being tested 
 // FlowMachine  cmpTest (redeclare package Medium = MyAir, 
  //                      isComp = true, eff = 0.95) ;  

  // Pressure controlling Compressor
  FlowMachine  cmpTrim (redeclare package Medium = MyAir, 
                    isComp = true, eff = 0.95) ;  

  // Cooler after compression
  HeatX hx (redeclare package Medium = MyAir, cf = 1.0e-3, 
                       A_wf = 1,  h_wf = 150, 
                       w_m = 1, w_cp = 420, holdup = 10e-3);

  // Blow valve Tee junction
  portMixer TeeBlow (redeclare package Medium = MyAir, 
                         N=3, Adiabatic = true) ;

  // Cooler outlet Tee junction
  portMixer TeeHx (redeclare package Medium = MyAir, 
                         N=3, Adiabatic = true) ;

  // DownStream Volume connector Tee
  // Test Compressor outlet and dowStream Trim feed join here
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
     //connect (volUp.outlet, cmpTest.inlet) ;

     // Downstream Test Chamber
     //connect (cmpTest.outlet, TeeComp.port[1]) ;
     connect (vDown.outlet, TeeComp.port[1]) ;
     connect (volDown.inlet, TeeComp.port[2]) ;
    
     connect (volDown.outlet, vBack.inlet) ;
     connect (vBack.outlet, atm.port) ;
      
     vBack.spo = 90 ;
     vBlow.spo = 30 ;
     vUp.spo = 90 ;
     vDown.spo = 2 ;

     // Limit power input as pressure in volumes reach set values
     if initial() then
            cmpTrim.Ws = 0 ;
     else
           cmpTrim.Ws = max (0, 100 * (1.7 - volUp.p/1e5)/1.7) ;
     end if ;
     cmpTrim.outlet.p =  (maxPr  - (maxPr - 1)
            *(cmpTrim.inlet.m_flow/maxFlow)) * cmpTrim.inlet.p ; 

/*
     cmpTest.Ws =  max (0, 10 * (1.1 - volDown.p/1e5)/1.1) ;
     cmpTest.outlet.p =  (maxPr  - (maxPr - 1)
             *(cmpTest.inlet.m_flow/maxFlow)) * cmpTest.inlet.p ; 
*/

     hx.Q_ew = -10 ; // External to wall heat transfer
     volUp.Q_loss = 30 ;
     volDown.Q_loss = 30 ;

initial algorithm
    volUp.T := 300 ;  // Initial Temperature
    volUp.p := 1e5  ;  // Initial Pressure
    volUp.Xi := Air ;

    volDown.T := 300 ;  // Initial Temperature
    volDown.p := 1e5  ;  // Initial Pressure
    volDown.Xi := Air ;
    volDown.m := volDown.vol * 1.0 ;

    hx.Tf := 300 ;  // Heater fluid Temp
    hx.Tw := 300 ;  // Heater wall Temp
//    hx.medium.Xi := Air ;
    
     vBlow.po :=25 ;
     vUp.po := 25 ;
     vDown.po := 25 ;
     vBack.po := 25 ;

     TeeBlow.p := 1e5 ;
     TeeBlow.T := 300 ;
     TeeBlow.Xi := Air ;

     TeeComp.p := 1e5 ;
     TeeComp.T := 300 ;
     TeeComp.Xi := Air ;

     TeeHx.p := 1e5 ;
     TeeHx.T := 300 ;
     TeeHx.Xi := Air ;

// compute consistent initial masses and expose medium starts numerically

volUp.medium.T := 300 ;
volUp.medium.p := 1e5 ;
volUp.medium.Xi := MyAir.reference_X[1:MyAir.nXi] ;

volDown.medium.T := 300 ;
volDown.medium.p := 1e5 ;
volDown.medium.Xi := MyAir.reference_X[1:MyAir.nXi] ;

volUp.m   := volUp.medium.d * volUp.vol;
volDown.m := volDown.medium.d * volDown.vol;
TeeBlow.m  := TeeBlow.medium.d * TeeBlow.vol;
TeeHx.m    := TeeHx.medium.d * TeeHx.vol;
TeeComp.m  := TeeComp.medium.d * TeeComp.vol;

//     cmpTrim.hout := MyAir.specificEnthalpy_pTX(1e5,300,MyAir.reference_X);
   //  cmpTest.hout := MyAir.specificEnthalpy_pTX(1e5,300,MyAir.reference_X);
    
end plant;
