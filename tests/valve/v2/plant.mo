model plant
/*
  Author: Ravi Saripalli
      Testing portMixer .... Works very well
*/

  import ThermoS.Uops.Reservoir;
  import ThermoS.Uops.Tanks.portMixer;
  import ThermoS.Uops.Valves.Valve ;
  import ThermoS.Uops.Valves.Vchar ;
  import ThermoS.Types.*;
  package Gas = ThermoS.Media.MyGas;
  constant Real AirComp[2] = {0.767,0.233};

  Reservoir     src	(redeclare  package Medium = Gas,
                               p = 3e5,  T = 300,  Xi = AirComp); // Reservoir 1

  Reservoir     sink	(redeclare  package Medium = Gas,
                               p = 1e5,  T = 300,  Xi = AirComp); // Reservoir 1

//  Valve v1 (redeclare each package Medium = Gas, vchar = Vchar.Linear);
//  Valve v2 (redeclare each package Medium = Gas, vchar = Vchar.EquiPercent);

  Valve v1 (redeclare each package Medium = Gas, cv=0.004/sqrt(0.5e5));
  Valve v2 (redeclare each package Medium = Gas, cv=0.004/sqrt(0.5e5), vchar=Vchar.EquiPercent);

  portMixer Node(redeclare package Medium = Gas,  vol=10e-6, N=3, Tset=400) ;

equation

       v1.po =   40  ;
       v2.po =   40  ;
       connect (src.port, Node.port[1]) ;
       connect (v1.inlet, Node.port[2]) ;
       connect (v2.inlet, Node.port[3]) ;
       connect(v2.outlet, sink.port) ;
       connect(v1.outlet, sink.port);
   
end plant;
