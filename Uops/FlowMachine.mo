within ThermoS.Uops;
model FlowMachine

/* 
     Adiabatic Compression/Turbine unit
     using Compressible Media
     Author: Ravi Saripalli  14th Nov. 2025
*/

  replaceable package Medium = PartialMixtureMedium ;

  FluidPort 	inlet (redeclare package Medium = Medium)  ; 
  FluidPort 	outlet (redeclare package Medium = Medium)  ; 

  Medium.ThermodynamicState	state_in    ; // Fluid state at inlet port
  Medium.ThermodynamicState	state_out   ; // Fluid state at outlet port
  Medium.ThermodynamicState	state_is    ; // Fluid state at isentropic 
  Medium.Temperature            Tin         ; // Machine Inlet Temp
  Medium.Temperature            Tis         ; // Isentropic Temp
  Medium.Temperature            Tout        ; // Compressor Outlet Temp
  Medium.MassFraction           X [Medium.nS] ; 
  SpecificEntropy               s ;  // Entropy of inlet/outlet stream
  SpecificEnthalpy              hin, hout, his;
  Power 	                Ws ; // Power delivered by Machine 

  parameter Boolean   isComp (start =  true) ;  // if false it is Expander 
  parameter Fraction  eff (start = 0.95) ;      //Isentropic Efficiency
  parameter Real      prat (start = if (isComp) then 1 else 0.9,
                           min = if (isComp) then 1 else 0.1,
                           max = if (isComp) then 5 else 1 ) ;

  equation

  // Mass balance 
   inlet.m_flow + outlet.m_flow = 0  ;     // No accumulation of mass

     outlet.Xi_outflow = inStream (inlet.Xi_outflow) ;
     inlet.Xi_outflow = outlet.Xi_outflow ;

     hin =  inStream (inlet.h_outflow) ;
     inlet.h_outflow = hin ;  // for completeness (revflow)
     hout =  outlet.h_outflow ;

     X  = cat(1, inlet.Xi_outflow, {1 - sum(inlet.Xi_outflow)});

     // Get outlet pressure 
     outlet.p =  prat * inlet.p ;

     // Get inlet state, entropy 
     state_in = Medium.setState_phX (inlet.p, hin, X);
     s = Medium.specificEntropy (state_in) ;
     Tin = Medium.temperature_phX (inlet.p, hin, X) ;

    // Determine outlet state if it were isentropic (state_iso)
     state_is = Medium.setState_psX (outlet.p, s, X);
     his = Medium.specificEnthalpy_pTX (outlet.p, Tis, X) ; 

    // Determine outlet Enthalpy and state 
     hout = hin + (his - hin) / eff ;
     state_out =  Medium.setState_phX (outlet.p, hout, X);
     Tout = Medium.temperature_phX (outlet.p, hout, X) ; 

    // Compressor Power Relation (we change work sign so power is >0 always
     Ws = inlet.m_flow * (hout - hin) * (if (isComp) then 1 else -1) ;
     
end FlowMachine;
