within Rig;

model Comp

/* 
     Adiabatic Compressor
     Author: Ravi Saripalli  14th Nov. 2025
*/

  replaceable package Medium = PartialMixtureMedium ;

  // Uni Directional Flow from inlet to outlet
  FluidPort  inlet (redeclare package Medium = Medium, m_flow (min = 0)) ;
  FluidPort outlet (redeclare package Medium = Medium, m_flow (max = 0)) ;

  Medium.ThermodynamicState	state_in    ; // Fluid state at inlet port
  Medium.ThermodynamicState	state_is    ; // Fluid state at isentropic 
  Medium.Temperature            Tin         ; // Machine Inlet Temp
  Medium.Temperature            Tis         ; // Isentropic Temp
  Medium.Temperature            Tout        ; // Compressor Outlet Temp
  SpecificEnthalpy              hin, hout, his;
  SpecificEntropy               s ;    // Entropy of inlet/outlet stream
  Power 	                Ws ;   // Power delivered by Machine 
  Medium.BaseProperties         medium ; // (preferredMediumStates = true) ;

  parameter Fraction  eff = 0.95     ;  //Isentropic Efficiency
  parameter Volume    vol = 10e-3    ; // 10 liters of holdup

  equation

   Tout = medium.T ;
   medium.Xi = inlet.Xi_outflow; 
   medium.p = outlet.p ;  // This is what I missed all along
   hout = medium.h ;
     
  // Mass balance 
   inlet.m_flow + outlet.m_flow =  0  ;     // No accumulation of mass

     hin = inStream (inlet.h_outflow) ;
     hout =  outlet.h_outflow ;
     inlet.h_outflow = medium.h ;

   // No composition Change
     inlet.Xi_outflow = inStream (inlet.Xi_outflow) ;
     outlet.Xi_outflow = inlet.Xi_outflow ; 
   
     // Get inlet state, entropy 
     state_in = Medium.setState_phX (inlet.p, hin, medium.Xi);
     s = Medium.specificEntropy (state_in) ;
     hin = Medium.specificEnthalpy_pTX (inlet.p, Tin, medium.Xi) ;

    // Determine outlet state if it were isentropic (state_iso)
     state_is = Medium.setState_psX (outlet.p, s, medium.Xi);
     his = Medium.specificEnthalpy_pTX (outlet.p, Tis, medium.Xi) ; 

    // Determine outlet Enthalpy and state 
     hout = hin + (his - hin) / eff ;
     hout = Medium.specificEnthalpy_pTX (outlet.p, Tout, medium.Xi) ; 

     Ws = inlet.m_flow  * (hout - hin) ;
     outlet.Xi_outflow = medium.Xi ;
end Comp;
