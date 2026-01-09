model plant
/*
  Author: Ravi Saripalli
  	16th May 2023 
*/
  import ThermoS.Types.*;
  import ThermoS.Media.MyGas;
  import ThermoS.Uops.Misc.Feed ;
  import ThermoS.Uops.Valves.Valve;
  import ThermoS.Uops.Tanks.GasTank ;
  import ThermoS.Uops.Tanks.Reservoir ;
  import ThermoS.Uops.Control.PIDControl ;

  constant    Real Air[MyGas.nXi] = {0.79, 0.21} ;
  Feed        supply (redeclare package Medium = MyGas); // InletFlow to tank
  Valve       valve (redeclare package Medium = MyGas , cv = (150e-3/60) / sqrt(1e3)) ;
  GasTank     tank (redeclare package Medium = MyGas, vol = 0.2 , Q_loss = 0); 
  Reservoir   atm (redeclare package Medium = MyGas, p=1e5, T=300, Xi=Air); // Reservoir 1
  PIDControl  pid (Kc = 0.5 , Ti = 10, Td = 0, reverseActing = true,
                   pvMin = 1e5, pvMax = 11e5, 
                   mvMin = 0, mvMax = 60e-3/60);
  
equation
     connect (supply.outlet, tank.inlet) ;
     connect (tank.outlet, valve.inlet) ;
     connect (valve.outlet, atm.port) ;

     pid.sp = 1.2e5 ; //1e5 + 6e5 * (1 - exp(-time/10)) ; // Tank Pressure setpoint
     pid.pv = tank.p ;
     pid.mv =  tank.inlet.m_flow ; 

     // supply.mdot = 1000 * 1e-3/60 ; // + 4.95 * sin(6*time) ; 
     supply.T = 300 ; // for a force feed you need flow, temp and comp
     supply.Xi = fill(1.0/MyGas.nS, MyGas.nXi) ;

     valve.po  = 10 * abs(sin(6*0.01*time))  ;
     
initial algorithm
    tank.medium.T := 300 ;  // Initial Temperature
    tank.medium.p := 1e5 ;  // Initial Temperature
    tank.medium.Xi := Air ;
end plant;
