model plant
/*
  Author: Ravi Saripalli
  	 9th May 2014 
*/
  import ThermoS.Uops.*;
  import ThermoS.Types.*;
  import Gas = ThermoS.Media.MyGas;

  constant Real AirComp[2] = {0.767,0.233}; // Note: reduced Xi is true

  Reservoir     res1	(redeclare package Medium = Gas, 
                         p = 5e5, T = 300, Xi = AirComp); // Reservoir 1
  Reservoir 	res3	(redeclare package Medium = Gas, 
                         p = 1e5, T = 300, Xi = AirComp); // Reservoir 2

  HeaterCooler htr(redeclare package Medium = Gas, cf = 1.0e-3, 
                                        A_wf = 1,  h_wf = 150, 
                                        w_m = 1, w_cp = 420, holdup = 50);


equation

    connect (res1.port, htr.inlet) ;
    connect (res3.port, htr.outlet) ;
    htr.Q_ew = 5e4 ;

initial algorithm
    htr.Tf :=300 ;
    htr.Tw :=300 ; 
end plant;
