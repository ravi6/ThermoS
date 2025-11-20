within ThermoS.Uops;
model FlowMachine

/* 
     Adiabatic Compression/Turbine unit
     using Compressible Media
     Author: Ravi Saripalli  14th Nov. 2025
*/

  replaceable package Medium = PartialMixtureMedium ;

  // Uni Directional Flow from inlet to outlet
  FluidPort  inlet (redeclare package Medium = Medium, m_flow (max = 0)) ;
  FluidPort outlet (redeclare package Medium = Medium, m_flow (min = 0)) ;

  Medium.ThermodynamicState	state_in    ; // Fluid state at inlet port
  Medium.ThermodynamicState	state_out   ; // Fluid state at outlet port
  Medium.ThermodynamicState	state_is    ; // Fluid state at isentropic 
  Medium.Temperature            Tin         ; // Machine Inlet Temp
  Medium.Temperature            Tis         ; // Isentropic Temp
  Medium.Temperature            Tout        ; // Compressor Outlet Temp
  SpecificEnthalpy              hin, hout, his;
  SpecificEntropy               s ;    // Entropy of inlet/outlet stream
  Power 	                Ws ;   // Power delivered by Machine 

  parameter Boolean   isComp (start =  true) ;  // if false it is Expander 
  parameter Fraction  eff (start = 0.95) ;      //Isentropic Efficiency

  equation

  // Mass balance 
   inlet.m_flow + outlet.m_flow = 0  ;     // No accumulation of mass

     hin = inStream (inlet.h_outflow) ;
     hout =  outlet.h_outflow ;
     outlet.Xi_outflow = Medium.reference_X[1:Medium.nXi] ;

     // Get inlet state, entropy 
     state_in = Medium.setState_phX (inlet.p, hin, Medium.reference_X);
     s = Medium.specificEntropy (state_in) ;
     Tin = Medium.temperature_phX (inlet.p, hin, Medium.reference_X) ;

    // Determine outlet state if it were isentropic (state_iso)
     state_is = Medium.setState_psX (outlet.p, s, Medium.reference_X);
     his = Medium.specificEnthalpy_pTX (outlet.p, Tis, Medium.reference_X) ; 

    // Determine outlet Enthalpy and state 
     hout = hin + (his - hin) / eff ;
     state_out =  Medium.setState_phX (outlet.p, hout, Medium.reference_X);
     Tout = Medium.temperature_phX (outlet.p, hout, Medium.reference_X) ; 

    // Compressor Power Relation (we change work sign so power is >0 always
    // The small derivative is to avoid algebraic loops
       0 = inlet.m_flow * (hout - hin) 
                           - Ws * (if (isComp) then 1 else -1) ;
     
end FlowMachine;
