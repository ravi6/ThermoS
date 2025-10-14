model plant
/*
  Author: Ravi Saripalli
  	11st Oct 2025 
    Compressor Tests 
*/
import ThermoS.Types.* ;
import ThermoS.Media.MyGas ;
import ThermoS.Uops.CompressorBasic ;
constant    Real Air[MyGas.nXi] = {0.7, 0.2} ;
CompressorBasic  comp (redeclare package Medium = MyGas, 
                                pr = 2,  n=1.4) ;
  Real x ;
equation
     comp.inlet.p = 1e5  ;
     comp.inlet.h_outflow = MyGas.specificEnthalpy( MyGas.setState_pTX(1e5, 300, Air));
    comp.inlet.Xi_outflow = Air ;
    comp.inlet.m_flow = 1;
    der(x) = 1 ; 

initial algorithm
     x := 10 ;
end plant;
