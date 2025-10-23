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
  SpecificEntropy               s ;   // Entropy of inlet/outlet stream
  Power 	                Ws ;  // Shaft work  (-ve for compressor)
  Energy                        U ;   // Holdup Internal Energy

  parameter Volume    holdup (start = 0.001) ;   // Compressor holdup
  parameter Fraction  eff (start = 1) ;          //Isentropic Efficiency

  equation

    // Mass balance 
     inlet.m_flow + outlet.m_flow = 0  ; // no mass accumulation

    // Energy balance 
     U = holdup * Medium.density(state_out) 
           * Medium.specificInternalEnergy(state_out) ;
     der(U) = - Ws + inlet.m_flow  * (hin - hout) ;

    // No change composition and flow is from inlet to outlet
     inlet.Xi_outflow = inStream(inlet.Xi_outflow) ; 
     outlet.Xi_outflow = inlet.Xi_outflow ;

     hin = inStream (inlet.h_outflow) ;
     inlet.h_outflow = hin ;  // A dummy assignment even though never used
     outlet.h_outflow = hout ;
     

     // Get inlet state, entropy 
     state_in = Medium.setState_phX (inlet.p, hin, inlet.Xi_outflow);
     s = Medium.specificEntropy (state_in) ;
     Tin = state_in.T ;

    // Determine outlet state if it were isentropic (state_iso)
     state_is = Medium.setState_psX (outlet.p, s, outlet.Xi_outflow);
     his = Medium.specificEnthalpy (state_is) ;

    // Determine outlet Enthalpy and state 
      hout = hin + (his - hin) / eff ;
      state_out = Medium.setState_phX (outlet.p, hout,
                          outlet.Xi_outflow);
      Tout = state_out.T ;

end CompressorBasic;
