within ThermoS.Uops;
model CompressorBasic

/* 
     Adiabatic Compression unit
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
  Medium.MassFraction           X_in  [Medium.nS], 
                                X_out [Medium.nS] ; // Full Composition
  SpecificEntropy               s ;  // Entropy of inlet/outlet stream
  SpecificEnthalpy              hin, hout, his;
  Power 	                Ws ; // Power delivered to comp.

  parameter Fraction  eff (start = 0.95) ;   //Isentropic Efficiency

  equation

    // Mass balance 
     inlet.m_flow + outlet.m_flow = 0  ;     // No accumulation of mass

     outlet.Xi_outflow = inStream (inlet.Xi_outflow) ;
     inlet.Xi_outflow = inStream (outlet.Xi_outflow) ;

     hin =  inStream (inlet.h_outflow) ;
     inlet.h_outflow = hin ;  // for completeness (revflow)
     hout =  outlet.h_outflow ;

     // stupid phX, psX state calls want full vector (so inconsistent)
     X_in  = cat(1, inlet.Xi_outflow, {1 - sum(inlet.Xi_outflow)});
     X_out = cat(1, outlet.Xi_outflow, {1 - sum(outlet.Xi_outflow)});
 

     // Get inlet state, entropy 
     state_in = noEvent (Medium.setState_phX (inlet.p, hin, X_in));
     s = noEvent (Medium.specificEntropy (state_in)) ;
     Tin = state_in.T ;

    // Determine outlet state if it were isentropic (state_iso)
     state_is = noEvent (Medium.setState_psX (outlet.p, s, X_out));
     his = noEvent (Medium.specificEnthalpy (state_is)) ;

    // Determine outlet Enthalpy and state 
     hout = hin + (his - hin) / eff ;
     state_out = noEvent (Medium.setState_phX (outlet.p, hout, X_out));
     Tout = state_out.T ;

   // Shaft work in adiabatic compression
     Ws = inlet.m_flow * (hout - hin) ;

end CompressorBasic;
