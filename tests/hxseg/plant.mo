model plant
/*
  Author: Ravi Saripalli
  	 9th May 2014 
*/
  import ThermoS.Uops.*;
  import ThermoS.Types.*;
  import ThermoS.Media.MyGas;

  constant Real AirComp[3] = {0.767,0.233,0};

  Reservoir     res1	(redeclare package Medium = MyGas, p = 5e5, T = 300, Xi = AirComp); // Reservoir 1
  Reservoir 	res2	(redeclare package Medium = MyGas, p = 5e5, T = 600, Xi = AirComp); // Reservoir 2
  Reservoir 	res3	(redeclare package Medium = MyGas, p = 1e5, T = 300, Xi = AirComp); // Reservoir 2
  Reservoir 	res4	(redeclare package Medium = MyGas, p = 1e5, T = 300, Xi = AirComp); // Reservoir 2
  HxSeg 	hx	(redeclare package Medium_c = MyGas, redeclare package Medium_h = MyGas); 

  HeaterCooler htr(redeclare package Medium = MyGas, cf = 5.0e-3, 
                                        A_wf = 1,  h_wf = 150, 
                                        w_m = 1, w_cp = 420, holdup = 50);

equation

     connect (res1.port, hx.portA_c);
     connect (res3.port, hx.portB_c);

     connect (res2.port, hx.portA_h);
     connect (res4.port, hx.portB_h);

    connect (res1.port, htr.inlet) ;
    connect (res3.port, htr.outlet) ;
    htr.Q_ew = 5e4 ;

initial algorithm
    hx.seg_h.Tf := 300 ;
    hx.seg_c.Tf := 300 ;
    hx.Tw := 300 ;

    htr.Tf :=300 ;
    htr.Tw :=300 ; 

end plant;
