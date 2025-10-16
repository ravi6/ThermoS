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
  Medium.ThermodynamicState	state_iso   ; // Fluid state at isentropic 
  Medium.Temperature            Tin (start = 300)        ; // Compressor Inlet Temp
  Medium.Temperature            Tout (start = 300)       ; // Compressor Outlet Temp

  HeatFlowRate                  Q (start=0) ; // Heat input to the device
  Power 	                Ws (start=0); // Power delivered to comp.

  Real 	effPoly (min = 0, max = 1, start = 1) ;   //Polytropic Efficiency
  Real 	effIs (min = 0, max = 1, start = 1);      //Isentropic  Efficiency
  Real  gamma (min = 0.1, max = 2, start = 1.4);  // Cp/Cv at inlet conditions

  parameter Real  n (min = 1, max = 1.5);              // Polytropic coeff
  parameter Real  pr (min = 0.1, max = 10, start = 1); // Pressure ratio 
  parameter Real  mdot (min = 0.0, max = 10,
                                     start = 0); // Pressure ratio 


  equation
    // Mass balance 
     inlet.m_flow + outlet.m_flow = 0  ;     // No accumulation of mass
     der (inlet.m_flow) = mdot ; 
     inStream(inlet.Xi_outflow) = outlet.Xi_outflow ;  // No change in gas comp 
     inStream(outlet.Xi_outflow) = inlet.Xi_outflow ;  // No change in gas comp 

     outlet.p = pr *  inlet.p ;
     gamma = 1.4 ; // Medium.isentropicExponent (state_in);
     effPoly = (1 - 1 / gamma)  / (1 - 1 / n) ;
     effIs =   (pr ^ ((1 - 1 / gamma) - 1)) 
               / (pr ^ ((1 - 1 / gamma) /  effPoly) - 1 ) ;

/*
    // Inlet state
     state_in = Medium.setState_phX (inlet.p, inStream(inlet.h_outflow),
                           inStream(inlet.Xi_outflow));
     Tin = state_in.T ;

    // Determine outlet state if it were isentropic (state_iso)
     state_is = Medium.setState_psX (outlet.p, Medium.specificEntropy (state_in),
                           inStream(outlet.Xi_outflow));

    // Exit Temperature when polyTropic compression
     Tout = state_in.T +  (state_is.T - state_in.T) / effPoly ;

    // Now  we have full state condition at outlet
     state_out = Medium.setState_pTX (outlet.p, Tout,
                          inStream(outlet.Xi_outflow));

   // Shaft work in polytropic compression
     Ws = inlet.m_flow * (Medium.specificEnthalpy (state_is) -
                          Medium.specificEnthalpy (state_in))
                         / effPoly ;

   //  Q Heat loss from the  compressor
   //  Ws Shaft work done on the compressor
   //  dH =  - Q +  Ws  (flow energy balance)
    Ws - Q   =  inlet.m_flow * (actualStream (inlet.h_outflow) +
                                actualStream (outlet.h_outflow)) ;
 */    
end CompressorBasic;
