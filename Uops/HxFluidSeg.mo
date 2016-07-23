within ThermoS.Uops;
model HxFluidSeg

  /* A generic Heat Exchanger Fluid Volume
     . No composition change 
     . Flow and pressure drop follow simple valve eqn.	
*/

  replaceable package Medium = PartialMixtureMedium ;
  FluidPort 	portA (redeclare package Medium = Medium)  ; 
  FluidPort 	portB (redeclare package Medium = Medium)  ; 

  parameter Real	cf	= 1     ; // Pressure Loss Coeff.      
  parameter Volume	holdup  = 1  ; // Heater fluid holdup 

  Medium.ThermodynamicState	state 		    ; // Fluid state (at actual outflow)
  Heat			        Qin	    	    ; // Heat input to the fluid
  Temperature			Tf		    ; // Heater Fluid temperature (K)
  Energy			U		    ; // Internal energy of fluid holdup
  SpecificHeatCapacity		Cp		    ; // Fluid Heat Capacity

  equation
    // Mass balance 
     portA.m_flow + portB.m_flow = 0  ;

    // Ignoring Composition change dynamics due to hold up
     portB.Xi_outflow = inStream(portA.Xi_outflow) ;  // Normal flow
     portA.Xi_outflow = inStream(portB.Xi_outflow) ;  // for  Reverse flow
     portA.h_outflow = portB.h_outflow ;   // Well mixed condition 

    /* Pressure differential drives flow 
       See valve eqn. for explanation  */
       portA.m_flow =   cf *  sqrt(Medium.density(state))
                           * regRoot((portA.p - portB.p), 10) ;
       state = Medium.setState_pTX((portA.p + portB.p)*0.5 , Tf, portB.Xi_outflow);
       portB.h_outflow = Medium.specificEnthalpy(state);

    // Energy Balance
       portA.m_flow * actualStream(portA.h_outflow) 
             + portB.m_flow *  actualStream(portB.h_outflow) + Qin 
        =  der(U)  ;   // fluid
        U = holdup * Medium.density(state) * Medium.specificInternalEnergy(state) ;
        Cp = Medium.specificHeatCapacityCp(state) ;
end HxFluidSeg;
