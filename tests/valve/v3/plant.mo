model plant
/*
  Author: Ravi Saripalli
*/

  import ThermoS.Uops.*;
  import ThermoS.Uops.Valves.*;
  import ThermoS.Types.*;
  package Gas = ThermoS.Media.MyGas(MassFlowRate(min=-10,max=20));
  constant Real AirComp[2] = {0.767,0.233};

  Reservoir     src	(redeclare  package Medium = Gas,
                               p = 2.0e5, each T = 300,  Xi = AirComp); 

  Reservoir     sink	(redeclare each package Medium = Gas,
                               p = 2.0e5,  T = 300,  Xi = AirComp); 

  Valve v1 (redeclare  package Medium = Gas,  cv=0.004/sqrt(0.5e5), Compressible = true) ;
  Valve v2 (redeclare  package Medium = Gas,  cv=0.004/sqrt(0.5e5), Compressible = true) ;


equation

       v1.po =    100  ;
       v2.po =    0   ;
       connect (src.port, v1.inlet) ;
       connect (src.port, v2.inlet) ;
       connect (v1.outlet, sink.port) ;
       connect (v2.outlet, sink.port) ;
   
end plant;
