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
  Medium.ThermodynamicState	state_out   ; // Fluid state at isentropic 
  Medium.Temperature            Tin         ; // Machine Inlet Temp
  Medium.Temperature            Tis         ; // Isentropic Temp
  Medium.Temperature            Tout        ; // Compressor Outlet Temp
  SpecificEnthalpy              hin, hout, his;
  SpecificEntropy               s ;    // Entropy of inlet/outlet stream
  Power 	                Ws ;   // Power delivered by Machine 

  parameter Fraction  eff = 0.95     ;  //Isentropic Efficiency
  parameter Real      pr  = 1.1      ;  //Pressure Ratio
  
  equation

  // Mass balance 
   inlet.m_flow + outlet.m_flow =  0  ;     // No accumulation of mass

   hin = inStream (inlet.h_outflow) ;
   hout =  outlet.h_outflow ;
   hin = inlet.h_outflow ;

    // Composition does not change
    outlet.Xi_outflow = inlet.Xi_outflow ; 
    inlet.Xi_outflow = Medium.reference_X [1:Medium.nXi] ;
    
   
     // Get inlet state, entropy 
     state_in = Medium.setState_phX (inlet.p, hin, inlet.Xi_outflow);
     s = Medium.specificEntropy (state_in) ;
     Tin = state_in.T ;
     //hin = Medium.specificEnthalpy_pTX (inlet.p, Tin, inlet.Xi_outflow) ;

    // Determine outlet state if it were isentropic (state_iso)
     state_is = Medium.setState_psX (outlet.p, s, outlet.Xi_outflow);
     Tis = state_is.T ;
     his = Medium.specificEnthalpy_pTX (outlet.p, Tis, outlet.Xi_outflow) ; 

    // Determine outlet Enthalpy and state 
     hout = hin + (his - hin) / eff ;
     state_out = Medium.setState_phX (outlet.p, hout, outlet.Xi_outflow);
     Tout = state_out.T ;

     Ws = inlet.m_flow  * (hout - hin) ;
     outlet.p = inlet.p * pr ;
end Comp;
