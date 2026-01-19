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
  SpecificEntropy               s  ;   // Entropy of inlet/outlet stream
  Power 	                Ws ;   // Power delivered by Machine 
  Real                          pr ;   // Pressure Ratio (output)

   // Hold up properties in Compressor
    Medium.BaseProperties       medium ;
    Enthalpy                    H      ; 

  parameter Volume    vol = 1e-4  ;     // holdup
  parameter Fraction  eff = 0.95     ;  //Isentropic Efficiency
  
  equation

   // Control volume state is outlet state
      medium.Xi = Medium.reference_X[1:Medium.nXi] ;
      medium.T = Tout ;  
      medium.p = outlet.p ;
      pr = outlet.p / inlet.p  ; 

  // Mass balance 
      inlet.m_flow + outlet.m_flow = 0   ;     // No accumulation of mass

   //  Energy balance
      H = vol * medium.d  * medium.h ;
      der(H)  = inlet.m_flow * (hin -  hout) + Ws/eff  + vol * der (medium.p) ;

   //  Well mixed and no change in composition
       inlet.Xi_outflow = medium.Xi ;
       outlet.Xi_outflow = medium.Xi ;
       outlet.h_outflow = hout ;
       inlet.h_outflow = hout ; 

   // Uni Directional model with  Control Volume
       hin = inStream (inlet.h_outflow) ;
   
   // Get inlet state, entropy 
       state_in = Medium.setState_phX (inlet.p, hin, medium.Xi);
       s = Medium.specificEntropy (state_in) ;
       Tin = state_in.T ;

    // Determine outlet state if it were isentropic (state_iso)
       state_is = Medium.setState_psX (outlet.p, s, medium.Xi);
       Tis = state_is.T  ;
       his = Medium.specificEnthalpy (state_is) ; 

    // Determine outlet Enthalpy and state 
       state_out = Medium.setState_phX (outlet.p, hout, medium.Xi);
       Tout = state_out.T ;

end Comp;
