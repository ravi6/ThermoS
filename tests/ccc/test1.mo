model plant
/*
    Compressor Tests 
 atm -> comp -> htr -> tank -> vout -> atm
*/

  import Modelica.Units.SI.MassFlowRate;
  import Rig.* ;

  constant    MassFlowRate  maxFlow = 200 * 1e-3 / 60 ;
  constant    Real CV = 200 * (1e-3/60) / sqrt (10000) ;

  Res        atm (redeclare package Medium = Air, p = 3e5, T = 300);
  Res        atm2 (redeclare package Medium = Air, p = 2e5, T = 300);
  Tank       tank (redeclare package Medium = Air, vol = 2, cf = 1e-4); 
  HeatX      htr (redeclare package Medium = Air, A_wf = 1,  h_wf = 150, 
                             cf = 1e-4, holdup = 5, w_m = 1, w_cp = 420);

//  Valve      vout (redeclare package Medium = Air, cv = CV);
//  Comp       comp (redeclare package Medium = Air, eff = 1) ;

equation
     connect (atm.port, tank.inlet) ;
     connect (tank.outlet, htr.inlet) ;
     connect (htr.outlet, atm2.port) ;
   // Assuming 2m2 surface area, 15 outside heattransfer coeff
     tank.Q_loss = 0 ; //15 * 2 * (tank.T - (15 + 273)) ;
     htr.Q_ew = 500 ;

initial algorithm
   
    tank.T := 300 ;  // Initial Temperature
    tank.p := 1e5  ;  // Initial Pressure
/*
    tank.medium.T := tank.T ;
    tank.medium.p := tank.p ;
    tank.medium.Xi := Air.reference_X; 
    tank.M := tank.vol * tank.medium.d ;
    tank.U := tank.M * tank.medium.u ;
*/
    htr.Tf := 1300 ;
    htr.Tw := 300 ; 
    htr.p := 1e5 ;
/*
    htr.medium.T := htr.Tf ;
    htr.medium.p := htr.p ;
    htr.medium.Xi := Air.reference_X; 
    htr.M := htr.holdup * htr.medium.d ;
    htr.H := htr.M * htr.medium.h ;
*/
end plant;
