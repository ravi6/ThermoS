model plant 
/*
  Test the effect of reversing flows with reservoirs
   connected to the bed. Localized Reservoir model with ability
   to dynamically change its pressure, compostion etc.

*/

  import Modelica.Media.Interfaces.Types.* ;
  import ThermoS.Uops.Valves.RealValve;
  import ThermoS.Uops.Reservoir;
  import ThermoS.Uops.Tanks.Cvolume;
  package Medium = ThermoS.Media.MyAir ;

 // package Medium = ThermoS.Media.MyAir ;
  constant MassFraction Xstart[1] = {0.79};
  Reservoir    src   (redeclare package Medium = Medium, 
                          p = 5e5, T = 300, Xi = Xstart); // Reservoir 1
  Reservoir    sink   (redeclare package Medium = Medium, 
                         p = 1e5, T = 300, Xi = Xstart); // Reservoir 2

  RealValve    v1 (redeclare package Medium = Medium, cv = 1/sqrt(0.5e5));
  RealValve    v2 (redeclare package Medium = Medium, cv = 1/sqrt(0.5e5));
  
  // By default the buffer is IsoThermal at T=300 (look at Cvolume.mo)
  Cvolume    buffer(redeclare package Medium = Medium, vol = 1);
equation

    connect (src.port, v1.inlet);
    connect (v1.outlet, v2.inlet);
    connect (v2.outlet, sink.port);
    connect (buffer.port, v2.inlet);

    v1.spo = 50.0 ;
    v2.spo = 30.0 ;
   
initial equation
     buffer.p = 1e5 ;
     buffer.Xi = Xstart;
     v1.po = 0 ;
     v2.po = 0 ;

end plant;
