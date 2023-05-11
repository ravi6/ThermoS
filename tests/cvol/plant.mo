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
  package Medium = ThermoS.Media.MyAir(Density(nominal=1.0), 
                                       MassFlowRate(nominal=1))  ;

  constant MassFraction Xstart[1] = {0.79};
  Reservoir    src   (redeclare package Medium = Medium, 
                          p = 5e5, T = 300, Xi = Xstart); // Reservoir 1
  Reservoir    sink   (redeclare package Medium = Medium, 
                         p = 1e5, T = 300, Xi = Xstart, port.m_flow(start=0.0)); // Reservoir 2

  RealValve    v1 (redeclare package Medium = Medium, cv=1/sqrt(0.5e5));
  RealValve    v2 (redeclare package Medium = Medium, cv=1/sqrt(0.5e5));
  

  Cvolume    buffer(redeclare package Medium 
              = Medium, port.m_flow(start=1.0, nominal=1e-3), vol=1);
equation

    connect (src.port, v1.inlet);
    connect (v1.outlet, v2.inlet);
    connect (v2.outlet, sink.port);
    connect (buffer.port, v2.inlet);

    v1.spo = 50.0 ;
    v2.spo = 30.0 ;
   
initial equation
     buffer.p = 1e5 ;
     buffer.T = 300 ;
     buffer.Xi = Xstart;
     v1.po = 0 ;
     v2.po = 0 ;

end plant;
