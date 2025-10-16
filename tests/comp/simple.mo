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
                              pr = 2,  n=1.4, mdot=1) ;

equation
     comp.inlet.p = 1e5 * time  ;
     comp.inlet.h_outflow = MyGas.specificEnthalpy( 
                       MyGas.setState_pTX(1e5, 300, Air));
    comp.inlet.Xi_outflow = Air ;
 
end plant;
