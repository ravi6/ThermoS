within ThermoS.Uops;
model CompressorBasic

/* 
     Polytropic Compression unit
     using Compressible Media
*/

  replaceable package Medium = PartialMixtureMedium ;

  FluidPort 	inlet (redeclare package Medium = Medium)  ; 
  FluidPort 	outlet (redeclare package Medium = Medium)  ; 

  Medium.ThermodynamicState	state_in    ; // Fluid state at inlet port
  Medium.ThermodynamicState	state_out   ; // Fluid state at outlet port
  Medium.ThermodynamicState	state_is    ; // Fluid state at isentropic 
  Medium.Temperature            Tin         ; // Compressor Inlet Temp
  Medium.Temperature            Tout        ; // Compressor Outlet Temp

  SpecificEnthalpy              hin, hout, his;
  Density                       rho_out  ;
  SpecificEntropy               s ;   // Entropy of inlet/outlet stream
  Power 	                Ws ;  // Shaft work  (-ve for compressor)
  Energy                        U ;   // Holdup Internal Energy

  parameter Volume    holdup (start = 0.001) ;   // Compressor holdup
  parameter Fraction  eff (start = 1) ;          //Isentropic Efficiency

  equation

    // Mass balance 
     inlet.m_flow + outlet.m_flow = 0  ; // no mass accumulation

    // No change composition and flow is from inlet to outlet
     inlet.Xi_outflow = inStream(inlet.Xi_outflow) ; 
     outlet.Xi_outflow = inStream(inlet.Xi_outflow) ;
     inlet.h_outflow = outlet.h_outflow ; // well mixed

     hin = actualStream (inlet.h_outflow) ;
     hout = actualStream (outlet.h_outflow) ; 
     
    // Energy balance (avoiding state record)
     rho_out = Medium.density_phX (outlet.p, hout, outlet.Xi_outflow) ;
     U = holdup * rho_out *  (hout - outlet.p / rho_out) ;
     der(U) = - Ws + inlet.m_flow  * (hin - hout) ;

     // Get inlet state, entropy 
     Tin = Medium.temperature_phX (inlet.p, hin, inlet.Xi_outflow) ;
     s = Medium.specificEntropy (
              Medium.setState_pTX (inlet.p, Tin, inlet.Xi_outflow));

    // Determine outlet state if it were isentropic (state_iso)
     his = Medium.specificEnthalpy (
                 Medium.setState_psX (outlet.p, s, outlet.Xi_outflow));

    // Determine outlet Enthalpy and state 
     hout = hin + (his - hin) / eff ;
     Tout = Medium.temperature_phX (outlet.p, hout, outlet.Xi_outflow);

     state_out = Medium.setState_pTX (outlet.p, Tout,
                         outlet.Xi_outflow);
     state_in = Medium.setState_pTX (inlet.p, Tin, inlet.Xi_outflow);
     state_is = Medium.setState_psX (outlet.p, s, outlet.Xi_outflow);
end CompressorBasic;
