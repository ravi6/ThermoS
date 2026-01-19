model plant
/*
    Compressor Tests 
 atm -> comp -> tank -> tee -> htr -> atm
                            -> atm
*/

  import Modelica.Units.SI.MassFlowRate;
  import Rig.* ;

  constant    MassFlowRate  maxFlow = 200 * 1e-3 / 60 ;
  constant    Real CV = 200 * (1e-3/60) / sqrt (10000) ;

  Res        atm (redeclare package Medium = Air, p = 1e5, T = 300);
  Tank       tank (redeclare package Medium = Air, vol = 2, cf = 1e-4); 
  HeatX      htr (redeclare package Medium = Air, A_wf = 1,  h_wf = 150, 
                             cf = 1e-4, holdup = 5, w_m = 1, w_cp = 420);
  Comp       comp (redeclare package Medium = Air, eff = 1) ;
  Tee        tee (redeclare package Medium = Air) ;
                        
//  Valve      vout (redeclare package Medium = Air, cv = CV);

equation
     connect (atm.port, comp.inlet) ;
     connect (comp.outlet, tank.inlet) ;
     connect (tank.outlet, tee.port[1]) ;
     connect (tee.port[2], htr.inlet) ; 
     connect (htr.outlet, atm.port) ;
     connect (tee.port[3], atm.port) ;
   // Assuming 2m2 surface area, 15 outside heattransfer coeff
     tank.Q_loss = 0 ; //15 * 2 * (tank.T - (15 + 273)) ;
     htr.Q_ew = 500 ;
     comp.Ws = 100;
     comp.pr = 1.1 ;
   

initial algorithm
    tank.T := 300 ;  // Initial Temperature
    tank.p := 1e5  ;  // Initial Pressure
    htr.Tf := 300 ;
    htr.Tw := 300 ; 
    htr.p := 1e5 ;
    comp.medium.p := 1e5 ;
    comp.medium.T := 300 ;
end plant;
