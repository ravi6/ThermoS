model plant
/*
  Author: Ravi Saripalli
  	 9th May 2014 
*/
  import ThermoS.Uops.*;
  import ThermoS.Types.*;
  import ThermoS.Media.JP8;


  Reservoir     res1	(redeclare package Medium = JP8, p = 5e5, T = 300); // Reservoir 1
  Reservoir 	res3	(redeclare package Medium = JP8, p = 1e5, T = 300); // Reservoir 2

  HeaterCooler htr(redeclare package Medium = JP8, cf = 1.0e-3, 
                                        A_wf = 1,  h_wf = 150, 
                                        w_m = 1, w_cp = 420, holdup = 0.5);


equation

    connect (res1.port, htr.inlet) ;
    connect (res3.port, htr.outlet) ;
    htr.Q_ew = 5e5 ;

initial algorithm
    htr.Tf :=300 ;
    htr.Tw :=300 ; 
end plant;
