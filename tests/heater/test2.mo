model plant
/*
  Author: Ravi Saripalli
  	 9th May 2014 
*/
  import ThermoS.Uops.*;
  import ThermoS.Uops.Valves.*;
  import ThermoS.Uops.Tanks.*;
  import ThermoS.Types.*;
  import Gas = ThermoS.Media.MyGas;

  constant Real AirComp[2] = {0.767,0.233}; // Note: reduced Xi is true

  Reservoir     res1	(redeclare package Medium = Gas, 
                         p = 5e5, T = 300, Xi = AirComp); // Reservoir 1
  Reservoir 	res3	(redeclare package Medium = Gas, 
                         p = 1e5, T = 300, Xi = AirComp); // Reservoir 2

  Valve v1 (redeclare each package Medium = Gas, cv=4e-3/sqrt(0.5e5));
  Valve v2 (redeclare each package Medium = Gas, cv=4e-3/sqrt(0.5e5));
  portMixer Node(redeclare package Medium = Gas,  vol=10e-6, N=3, Adiabatic=true) ;
  HeaterCooler htr(redeclare package Medium = Gas, cf = 1.0e-3, 
                                        A_wf = 1,  h_wf = 150, 
                                        w_m = 1, w_cp = 420, holdup = 50);
equation

    connect (res1.port, v1.inlet) ;
    connect (res1.port, v2.inlet) ;
    connect (v1.outlet, Node.port[1]);
    connect (v2.outlet, Node.port[2]);
    connect (Node.port[3], htr.inlet);
    connect (htr.outlet, res3.port) ;

    v1.po = 50 ;
    v2.po = 80 ;
    htr.Q_ew = 5e4 ;

initial algorithm
    htr.Tf :=300 ;
    htr.Tw :=300 ; 
end plant;
