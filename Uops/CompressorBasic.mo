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
  Medium.Temperature            Tin (start = 300)        ; // Compressor Inlet Temp
  Medium.Temperature            Tout (start = 300)       ; // Compressor Outlet Temp

  SpecificEntropy               s ;  // Entropy of inlet/outlet stream
  SpecificEnthalpy              hin, hout, his;
  Power 	                Ws ; // Power delivered to comp.

  parameter Real  eff (min = 0, max = 1, start = 1) ;   //Isentropic Efficiency

  equation

    // Mass balance 
     inlet.m_flow + outlet.m_flow = 0  ;     // No accumulation of mass
     inStream(inlet.Xi_outflow) = outlet.Xi_outflow ;  
     outlet.Xi_outflow = inlet.Xi_outflow ;  
     hin =  (inlet.h_outflow) ;
     hout =  (outlet.h_outflow) ;

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

   // Shaft work in adiabatic compression
     Ws = inlet.m_flow * (hout - hin) ;

end CompressorBasic;
