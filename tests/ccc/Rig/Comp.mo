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
  Medium.ThermodynamicState	state_out   ; // Fluid state at isentropic 
  Medium.Temperature            Tin         ; // Machine Inlet Temp
  Medium.Temperature            Tout        ; // Compressor Outlet Temp
  SpecificEnthalpy              hin, hout, his;
  Power 	                Ws ;   // Power delivered by Machine 
  Real                          pr ;   // Pressure Ratio (output)

  parameter Fraction  eff = 0.95     ;  //Isentropic Efficiency
  
  equation

      pr = outlet.p / inlet.p  ; 

  // Mass balance 
      inlet.m_flow + outlet.m_flow = 0   ;     // No accumulation of mass

   //  Energy balance
       hout = hin + (his - hin) / eff ;
       Ws = inlet.m_flow * (hout - hin) ; 

   //  Well mixed and no change in composition
       outlet.Xi_outflow = Medium.reference_X ;
       inlet.Xi_outflow = outlet.Xi_outflow ;
       outlet.h_outflow = hout ;
       inlet.h_outflow = hout ; 

   // Uni Directional model with  Control Volume
       hin = inStream (inlet.h_outflow) ;
   
   // Get inlet state, entropy 
       state_in = Medium.setState_phX (inlet.p, hin, inlet.Xi_outflow);
       Tin = state_in.T ;

    // Determine outlet state if it were isentropic (state_iso)
       his = Medium.isentropicEnthalpy(outlet.p, state_in);

    // Determine outlet Enthalpy and state 
       state_out = Medium.setState_phX (outlet.p, hout, outlet.Xi_outflow);
       Tout = state_out.T ;

end Comp;
