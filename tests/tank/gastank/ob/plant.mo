model plant
/*
  Author: Ravi Saripalli
  	21st May 2023 
     Obogs Air Supply System Evaluation 
*/
  import ThermoS.Types.*;
  import ThermoS.Media.MyGas;
  import ThermoS.Uops.Feed ;
  import ThermoS.Uops.Valves.Valve;
  import ThermoS.Uops.Tanks.GasTank ;
  import ThermoS.Uops.Reservoir ;
  import ThermoS.Uops.OnOff ;
  import ThermoS.Uops.Controller ;

  constant    Real Air[MyGas.nXi] = {0.79, 0.21} ;

  Feed        supply (redeclare package Medium = MyGas); // InletFlow to tank
  Valve       valve (redeclare package Medium = MyGas , cv = (1500e-3/60) / sqrt(1e4)) ;
  GasTank     tank (redeclare package Medium = MyGas, vol = 0.2 , Q_in=0); 
  Reservoir   atm (redeclare package Medium = MyGas, p=1e5, T=300, Xi=Air); // Reservoir 1

 // Compressor Control
  OnOff   onoff (pvMin = 0, pvMax = 10e5, mvMin = 0, mvMax = 1000e-3/60,
                      deadBand = 0.5e5);

 // Obogs Demand Control
  Controller  pid ( pvMin = 0, pvMax = 1000e-3/60, 
                    mvMin = 0, mvMax = 100,
                    Kc =1 , Ti = 1, Td = 0,
                    reverseActing = true );
  
equation
     connect (supply.outlet, tank.inlet) ;
     connect (tank.outlet, valve.inlet) ;
     connect (valve.outlet, atm.port) ;

     onoff.sp = 1e5 + 6e5 * (1 - exp(-time/10)) ; // soft ramp up Tank Pressure setpoint
     onoff.pv = tank.p ;  // Controller Measure Value
     onoff.mv = tank.inlet.m_flow;

     // supply.mdot = 1000 * 1e-3/60 ; // + 4.95 * sin(6*time) ; 
     supply.T = 300 ; // for a force feed you need flow, temp and comp
     supply.Xi = fill(1.0/MyGas.nS, MyGas.nXi) ;

     //valve.po = 0 ;
     pid.mv = valve.po ; 
     pid.sp = (600e-3/60) * (1 - exp(-time)) ;
     pid.pv = valve.inlet.m_flow ;
     
initial algorithm
    tank.T := 300 ;  // Initial Temperature
    tank.p := 1e5 ;  // Initial Temperature
    tank.Xi := Air ;

end plant;
