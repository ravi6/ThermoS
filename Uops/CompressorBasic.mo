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
  Medium.ThermodynamicState	state_out   ;// Fluid state at outlet port
  Heat			        Q     	    ; // Heat input to the device
  Power 		        Ws     	    ; // Power delivered to the fluid by the shaft 

  Real 	effPoly(min = 0, max = 1) = 0.9;    //Polytropic Efficiency
  Real 	effIsen(min = 0, max = 1);          //Isentropic  Efficiency
  Real  n(min = 1, max = 1.5);              // Polytropic coeff
  Real  pr(min=0.1, max=10, start=1);       // Pressure ratio (Outlet/Inlet)
  Real  gamma (min=0.1, max=2, start=1.4);  // Cp/Cv at inlet conditions

  equation
    // Mass balance 
     inlet.m_flow + outlet.m_flow = 0  ;     // No accumulation of mass
     inlet.Xi_outflow = outlet.Xi_outflow ;  // No change in gas comp 

    // Isentropic Efficiency and Polytropic Efficiency relation
     state_in = Medium.setState_pTX (inlet.p, inlet.T, 
                           inStream(inlet.Xi_outflow));
     gamma = MoistAir.isentropicExponent(state_in);
     pr    = outlet.p / inlet.p ;
     effPoly = (gamm - 1) * n / (gamma * (n -1)) ;
     effIsen = ( (pr) ^ ((gamma - 1) / gamma) - 1 ) /
               ( (pr) ^ ((gamma - 1) / (gamma * effPoly)) - 1 ) ;

    // Outlet Temp after accounting for Isentropic efficincy
     outlet.T = pr ^ ((gamma - 1) / gamma / effIsen) ;
     state_out = Medium.setState_pTX (outlet.p, outlet.T, 
                                      outlet.Xi_outflow);
 
   // Work done is R(T2 - T1)/ (1-n) &  R = Cp - Cv 
     Ws  =    inlet.m_flow * ((gamma - 1)/(1 - n)) 
                         * (Medium.specificInternal (state_out) 
                          - Medium.specificEnthalpy (state_in)) ;
   //  Q = dH + Ws  (flow energy balance)
     Q  =  inlet.m_flow * (Medium.specificEnthalpy (state_out)
                          - Medium.specificEnthalpy (state_in))
           - Ws ;
     
end CompressorBasic;
