within ThermoS.Uops;
model CompressorBasic

  /* 
     Polytropic Compression unit
     using Compressible Media
     // Algebraic ... no derivatives :))
*/

  replaceable package Medium = PartialMixtureMedium ;

  FluidPort 	inlet (redeclare package Medium = Medium)  ; 
  FluidPort 	outlet (redeclare package Medium = Medium)  ; 

  Medium.ThermodynamicState	state_in    ; // Fluid state at inlet port
  Medium.ThermodynamicState	state_out   ; // Fluid state at outlet port
  HeatFlowRate                  dQ (start=0); // Heat input to the device
  Power 		        Ws (start=0); // Power delivered to comp.

  Real 	effPoly(min = 0, max = 1, start = 1) ;    //Polytropic Efficiency
  Real 	effIsen(min = 0, max = 1, start = 1);     //Isentropic  Efficiency
  Real  gamma (min=0.1, max=2, start=1.4);        // Cp/Cv at inlet conditions

  parameter Real  n(min = 1, max = 1.5);              // Polytropic coeff
  parameter Real  pr(min=0.1, max=10, start=1);       // Pressure ratio (Outlet/Inlet)

  equation
    // Mass balance 
     inlet.m_flow + outlet.m_flow = 0  ;     // No accumulation of mass
     inlet.Xi_outflow = outlet.Xi_outflow ;  // No change in gas comp 

    // Isentropic Efficiency and Polytropic Efficiency relation
     state_in = Medium.setState_phX (inlet.p, inStream(inlet.h_outflow),
                           inStream(inlet.Xi_outflow));
     state_out = Medium.setState_phX (outlet.p, inStream(outlet.h_outflow),
                                      inStream(outlet.Xi_outflow));

     gamma = Medium.isentropicExponent(state_in);
     outlet.p = pr *  inlet.p ;
     effPoly = (gamma - 1) * n / (gamma * (n -1)) ;
     effIsen = ( (pr) ^ ((gamma - 1) / gamma) - 1 ) /
               ( (pr) ^ ((gamma - 1) / (gamma * effPoly)) - 1 ) ;

    // Outlet Temp after accounting for Isentropic efficincy
     state_out.T = state_in.T * pr ^ ((gamma - 1) / gamma / effIsen) ;
 
   // Work done is R(T2 - T1)/ (1-n) &  R = Cp - Cv 
     Ws  =    inlet.m_flow * ((gamma - 1)/(1 - n)) 
                         * (Medium.specificInternalEnergy (state_out) 
                          - Medium.specificInternalEnergy (state_in)) ;
   //   dH = - dQ +  Ws  (flow energy balance)
   //  dQ > 0 Heat loss from compressor
   //  Ws > 0 Power input to the compressor
     dQ  =  inlet.m_flow * (Medium.specificEnthalpy (state_out)
                          - Medium.specificEnthalpy (state_in))
           - Ws ;
     
end CompressorBasic;
