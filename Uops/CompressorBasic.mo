within ThermoS.Uops;
model CompressorBasic

/* 
     Polytropic Compression unit
     using Compressible Media
*/

  replaceable package Medium = PartialMixtureMedium ;

  FluidPort 	inlet (redeclare package Medium = Medium)  ; 
  FluidPort 	outlet (redeclare package Medium = Medium)  ; 

  Medium.Temperature            Tin         ; // Compressor Inlet Temp
  Medium.Temperature            Tout        ; // Compressor Outlet Temp
  Medium.MassFraction           Xi[Medium.nXi]    ; // Medium composition

  SpecificEnthalpy              hin, hout, his;
  SpecificEntropy               s ;   // Entropy of inlet/outlet stream
  Power 	                Ws ;  // Shaft work  (-ve for compressor)
  
  parameter Fraction  eff (start = 1) ;          //Isentropic Efficiency

  equation

    // Mass balance 
     inlet.m_flow + outlet.m_flow = 0  ; // no mass accumulation

    // Energy balance 
     Ws = - inlet.m_flow  * (hout - hin) ;

     inlet.Xi_outflow = outlet.Xi_outflow ;
     outlet.Xi_outflow = inStream(inlet.Xi_outflow) ; 

     hin = actualStream (inlet.h_outflow) ;
     hout = actualStream (outlet.h_outflow) ;
     inStream (inlet.h_outflow) = hin ;
     
     // Get inlet Entropy 
     s = Medium.specificEntropy(
             Medium.setState_pTX(inlet.p, Tin, Xi));

    // Determine outlet isentropic Enthalpy
     his = Medium.specificEnthalpy_psX(outlet.p, s, Xi);

    // Determine outlet Enthalpy and state 
     hout = hin + (his - hin) / eff ;
     hin = Medium.specificEnthalpy  (Medium.setState_pTX(inlet.p, Tin, Xi));
     hout = Medium.specificEnthalpy (Medium.setState_pTX(outlet.p, Tout, Xi));

end CompressorBasic;
