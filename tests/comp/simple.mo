model plant
/*
  Author: Ravi Saripalli
  	11st Oct 2025 
    Compressor Tests 
*/
import ThermoS.Types.* ;
import ThermoS.Media.MyGas ;
import ThermoS.Uops.CompressorBasic ;
import ThermoS.Uops.Reservoir ;
import ThermoS.Uops.Valves.Valve ;
import ThermoS.Uops.Tanks.GasTank ;

constant    Real Air[MyGas.nXi] = {0.7, 0.2} ;
//CompressorBasic  comp (redeclare package Medium = MyGas, 
//                                          pr = 1.2,  n=1.4) ;
Reservoir   suction  (redeclare package Medium = MyGas, p=1.2e5, T=300, Xi=Air);
Reservoir   atm      (redeclare package Medium = MyGas, p=1e5, T=300, Xi=Air);
Valve       valve    (redeclare package Medium = MyGas , 
                                        cv = (1500e-3/60) / sqrt(4e5)) ;
GasTank     tank (redeclare package Medium = MyGas, vol = 0.2 , Q_in=0); 

equation
//    connect(suction.port, comp.inlet) ;
//    connect(comp.outlet, tank.inlet) ;
    connect(suction.port, tank.inlet) ;
    connect(tank.outlet, valve.inlet) ;
    connect(valve.outlet, atm.port);
   // comp.inlet.m_flow = 1.0 ;
    valve.po = 50 ;

initial algorithm
    tank.p := 1e5 ;
    tank.T := 300 ;
    tank.Xi := Air ;
end plant;
