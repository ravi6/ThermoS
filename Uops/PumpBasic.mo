within ThermoS.Uops;
model PumpBasic

  /* A simple Pump Model for incompressible fluids 

*/

  replaceable package Medium = PartialMixtureMedium ;

  FluidPort 	inlet (redeclare package Medium = Medium)  ; 
  FluidPort 	outlet (redeclare package Medium = Medium)  ; 

  parameter Volume		holdup  = 1e-3      ; // Pump fluid holdup 
  Medium.ThermodynamicState	state 		    ; // Fluid state in the pump
  Heat			        Q 	    	    ; // Heat input to the device
  Temperature			T		    ; // Heater Fluid temperature (K)
  Energy			U		    ; // Internal energy of fluid holdup

  Power 		Ws	         	    ; // Power delivered to the fluid by the shaft 
  Real 			eff(min = 0, max = 1) = 1   ; //(adiabatic efficiency)
  Length                head                        ; // 
  VolumeFlowRate        Vdot                        ;

  equation
    // Mass balance 
     inlet.m_flow + outlet.m_flow = 0  ;

    // Ignoring Composition change dynamics due to hold up
     outlet.Xi_outflow = inStream(inlet.Xi_outflow) ;  // Normal flow
     inlet.Xi_outflow = inStream(outlet.Xi_outflow) ;  // for  Reverse flow
     inlet.h_outflow = outlet.h_outflow ;   // Well mixed condition 

    /* State of the Control Volume */ 
       state = Medium.setState_pTX((inlet.p + outlet.p)*0.5 , T, outlet.Xi_outflow);
       outlet.h_outflow = Medium.specificEnthalpy(state);

    // Shaft work that is transferred to fluid
        head = (outlet.p - inlet.p) / (Medium.density(state) * 9.81) ;
        Ws =   Vdot * (outlet.p - inlet.p)  * eff ;
        inlet.m_flow = Vdot * Medium.density(state) ;

    // Energy Balance
       inlet.m_flow * actualStream(inlet.h_outflow)
             + outlet.m_flow *  actualStream(outlet.h_outflow) + Q - Ws 
        =  der(U)  ;   // fluid
        U = holdup * Medium.density(state) * Medium.specificInternalEnergy(state) ;
   
end PumpBasic;
