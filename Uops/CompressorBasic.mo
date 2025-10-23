within ThermoS.Uops;
model CompressorBasic

/* 
     Polytropic Compression unit
     using Compressible Media
*/

  replaceable package Medium = PartialMixtureMedium ;

  FluidPort 	inlet (redeclare package Medium = Medium)  ; 
  FluidPort 	outlet (redeclare package Medium = Medium)  ; 

  Medium.ThermodynamicState	state_out   ; // Fluid state at outlet port
  Medium.Temperature            Tin         ; // Compressor Inlet Temp
  Medium.Temperature            Tout        ; // Compressor Outlet Temp

  SpecificEnthalpy    hin, hout ; // inlet and outlet enthalpies
  Power 	      Ws ;  // Shaft work  (-ve for compressor)
  Energy              U ;   // Holdup Internal Energy

  Real    prat (min = 0.1, max = 10, start = 1) ; // pressure ratio

  parameter Volume    holdup (start = 0.001) ;   // Compressor holdup
  parameter Fraction  eff (start = 1) ;          //Isentropic Efficiency
  constant  Real      gamma (start = 1.4) ;

  equation

    // Mass balance 
     0  = inlet.m_flow + outlet.m_flow  ; // no mass accumulation

    // Energy balance 
     U = holdup * Medium.density(state_out) * Medium.specificInternalEnergy(state_out) ;
     der(U) = - Ws + inlet.m_flow  * (hin - hout) ;

    // Ignoring Composition change dynamics due to hold up
     outlet.Xi_outflow = inStream(inlet.Xi_outflow) ;  // Normal flow
     inlet.Xi_outflow = inStream(outlet.Xi_outflow) ;  // for  Reverse flow

     hin = inStream (inlet.h_outflow) ;
     Tin = Medium.temperature_phX(inlet.p, hin, inlet.Xi_outflow ) ;
     prat  = outlet.p / inlet.p  ;

     // Adiabatic 
     Tout = Tin * ( 1 +  (1/eff) * ( prat ^ (1 - 1/gamma)  - 1) ) ;
     hout = Medium.specificEnthalpy_pTX(outlet.p, Tout, outlet.Xi_outflow);
     hout =  outlet.h_outflow ;
     state_out = Medium.setState_pTX (outlet.p, Tout, outlet.Xi_outflow);
end CompressorBasic;
