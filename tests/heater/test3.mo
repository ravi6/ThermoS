model plant
/*
  Author: Ravi Saripalli
  	 9th May 2014 
  res1 -> v1 an v2  -> Node -> htr  -> v3 -> res3
                                    -> res3 
*/
  import ThermoS.Uops.*;
  import ThermoS.Uops.Valves.*;
  import ThermoS.Uops.Tanks.*;
  import ThermoS.Uops.HeatExch.*;
  import ThermoS.Types.*;
  import Gas = ThermoS.Media.MyGas;

  constant Real AirComp[2] = {0.767,0.233}; // Note: reduced Xi is true

  Reservoir     res1	(redeclare package Medium = Gas, 
                         p = 1.1e5, T = 300, Xi = AirComp); // Reservoir 1
  Reservoir 	res3	(redeclare package Medium = Gas, 
                         p = 1e5, T = 300, Xi = AirComp); // Reservoir 2

  Valve v1 (redeclare each package Medium = Gas, cv=4e-3/sqrt(0.5e2));
  Valve v2 (redeclare each package Medium = Gas, cv=4e-3/sqrt(0.5e2));
  portMixer Node(redeclare package Medium = Gas,  vol=1e-3, N=3, Adiabatic=true) ;
  HeatX htr(redeclare package Medium = Gas, cf = 1.0e-3, 
                                        A_wf = 1,  h_wf = 150, 
                                        w_m = 1, w_cp = 420, holdup = 1e-1);
  portMixer Node2 (redeclare package Medium = Gas,  vol=1e-3, N=3, Adiabatic=true) ;
  Valve v3 (redeclare each package Medium = Gas, cv=4e-3/sqrt(0.5e2));
//  Valve v4 (redeclare each package Medium = Gas, cv=4e-3/sqrt(0.5e2));
equation

    connect (res1.port, v1.inlet) ;
    connect (res1.port, v2.inlet) ;
    connect (v1.outlet, Node.port[1]);
    connect (v2.outlet, Node.port[2]);
    connect (Node.port[3], htr.inlet);
    connect (htr.outlet, Node2.port[1]);
    connect (v3.inlet, Node2.port[2]);
    connect (v3.outlet, res3.port);
    connect (Node2.port[3], res3.port);
    
   
    v1.po = 50 ;
    v2.po = 80 ;
    v3.po = 30 ;
    htr.Q_ew = 50000 ;

initial equation
    htr.Tf =300 ;
    htr.Tw =300 ; 
    htr.medium.p = 1e5 ;
  
end plant;
