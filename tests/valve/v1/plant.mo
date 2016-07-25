model plant
/*
  Author: Ravi Saripalli
*/

  import ThermoS.Uops.*;
  import ThermoS.Uops.Valves.*;
  import ThermoS.Types.*;
  package Gas = ThermoS.Media.MyGas;
  constant Real AirComp[2] = {0.767,0.233};

  Reservoir     src[1]	(redeclare each package Medium = Gas,
                              each p = 3e5, each T = 300, each Xi = AirComp); // Reservoir 1

  Reservoir     sink[1]	(redeclare each package Medium = Gas,
                              each p = 1e5, each T = 300, each Xi = AirComp); // Reservoir 1

  Valve valve[3] (redeclare each package Medium = Gas, 
                     //each inlet.Xi_outflow(start=AirComp),
                     //each outlet.Xi_outflow(start=AirComp),
                     //vchar = {Vchar.Linear, Vchar.Linear, Vchar.Linear}) ;
                     vchar = {Vchar.Linear, Vchar.EquiPercent, Vchar.FastActing}) ;

/* key word <each> is a short cut to repeated element array of the size of the model array
   being initialized but located before the variable ... make it read like english??*/


equation

    for k in 1:3 loop
       valve[k].po =   10 * time ;
       connect (src[1].port, valve[k].inlet) ;
       connect (valve[k].outlet, sink[1].port) ;
    end for;
   
end plant;
