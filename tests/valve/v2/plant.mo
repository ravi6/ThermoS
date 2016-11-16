model plant
/*
  Author: Ravi Saripalli
      Testing P3Mixer .... Works very well
*/

  import ThermoS.Uops.Reservoir;
  import ThermoS.Uops.Tanks.P3Mixer;
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

  Valve v1 (redeclare each package Medium = Gas, cv=0.004/sqrt(0.5e5), Compressible=true, dpTol=0.1);
  Valve v2 (redeclare each package Medium = Gas, cv=0.004/sqrt(0.5e5), vchar=Vchar.EquiPercent, Compressible=false);

  P3Mixer Node(redeclare package Medium = Gas,  vol=1.0e-6) ;

equation

       v1.po =   10 * time ;
       v2.po =   10 * time ;
       connect (src.port, Node.port1) ;
       connect (v1.inlet, Node.port2) ;
       connect (v2.inlet, Node.port3) ;
       connect(v2.outlet, sink.port) ;
       connect(v1.outlet, sink.port);
   
end plant;
