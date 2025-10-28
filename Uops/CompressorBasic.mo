within ThermoS.Uops;
model CompressorBasic

/* 
     Polytropic Compression unit
     using Compressible Media
     This version avoids calls to setState_phX, pTX, psX
     ... differentiability issues and derving T from H issues
*/

  replaceable package Medium = PartialMixtureMedium ;

  FluidPort 	inlet (redeclare package Medium = Medium)  ; 
  FluidPort 	outlet (redeclare package Medium = Medium)  ; 

  Medium.BaseProperties	        state_in    ; // Fluid state at inlet port
  Medium.BaseProperties	        state_out   ; // Fluid state at outlet port
  Medium.BaseProperties	        state_is    ; // Fluid state at isentropic 
  Medium.SpecificEntropy        s           ;
  Medium.Temperature            Tin         ; // Compressor Inlet Temp
  Medium.Temperature            Tout        ; // Compressor Outlet Temp

  SpecificEnthalpy              hin, hout, his;
  Power 	                Ws ;  // Shaft work  (-ve for compressor)
  

//  parameter Volume    holdup (start = 0.001) ;   // Compressor holdup
  parameter Fraction  eff (start = 1) ;          //Isentropic Efficiency

  equation

    // Mass balance (no holdup)
     inlet.m_flow + outlet.m_flow = 0 ;

    // We have no composition changes in the unit
     outlet.Xi_outflow = inlet.Xi_outflow ;  
     inStream (inlet.Xi_outflow) = outlet.Xi_outflow ;

    // Energy Balance 
     hin = actualStream (inlet.h_outflow) ;
     hout = actualStream (outlet.h_outflow) ; 
     Ws = - (inlet.m_flow * hin + outlet.m_flow * hout) ;

     // Get inlet state, entropy (avoid call to  setStat_phX) 
     state_in.Xi = inlet.Xi_outflow ;
     state_in.p = inlet.p ;
     state_in.h = hin ;
     s = Medium.specificEntropy (state_in.state) ;
     state_in.T = Tin ;

    // Determine outlet state if it were isentropic (state_iso)
     state_is.Xi = inlet.Xi_outflow ;
     state_is.p  = outlet.p ;
     state_is.h = Medium.specificEnthalpy (state_is.state) ;
     s = Medium.specificEntropy (state_is.state) ;

    // Determine outlet Enthalpy and state 
     hout = hin + (state_is.h - hin) / eff ;
     state_out.Xi = inlet.Xi_outflow ;
     state_out.p = outlet.p ;
     state_out.h = hout ; 
     state_out.T = Tout ;
end CompressorBasic;
