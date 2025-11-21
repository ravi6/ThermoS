model plant
/*
  Author: Ravi Saripalli
      Testing portMixer .... Works very well
*/

  import ThermoS.Types.*;
  import ThermoS.Uops.Tanks.Reservoir;
  import ThermoS.Uops.Tanks.portMixer;
  import ThermoS.Uops.Valves.RealValve ;
  package Gas = ThermoS.Media.MyGas;
  constant Real AirComp[2] = {0.767,0.233};
  constant Real CV = (200e-3/60)/sqrt(100);

  Reservoir     src	(redeclare  package Medium = Gas,
                               p = 1.30e5,  T = 300,  Xi = AirComp); // Reservoir 1
  Reservoir     sink	(redeclare  package Medium = Gas,
                               p = 1e5,  T = 300,  Xi = AirComp); // Reservoir 1

  RealValve v1 (redeclare  package Medium = Gas, cv=CV, tau=5);
  RealValve v2 (redeclare  package Medium = Gas, cv=CV, tau=5);
  RealValve v3 (redeclare  package Medium = Gas, cv=CV, tau=5);
  portMixer Node (redeclare package Medium = Gas,  vol=10e-6, N=3, Adiabatic=true) ;
  portMixer Node2 (redeclare package Medium = Gas,  vol=10e-6, N=3, Adiabatic=true) ;

equation

     v1.spo =  40  ;
     v2.spo =  40  ;
     v3.spo =  40  ;

     connect (src.port, Node.port[1]) ;
     connect (v1.inlet, Node.port[2]) ;
     connect (v2.inlet, Node.port[3]) ;
 
     connect (v1.outlet, Node2.port[1]) ;
     connect (v2.outlet, Node2.port[2]) ;

     connect(Node2.port[3], v3.inlet) ;
     connect(v3.outlet, sink.port) ;

initial algorithm
     v1.po := 10 ;
     v2.po := 10 ;   
     v3.po := 10 ;
end plant;
