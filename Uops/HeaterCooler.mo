// Version 1.3
//  24th Mar. 2014
within ThermoS.Uops;
model HeaterCooler

  /* This heater/cooler Model 
     . Fluid is completely mixed in the heater (CSTR)
     . Permits pressure loss specification through loss coeff
       based on volumetric flow, or as a fixed value.
     . Accounts for  thermal inertia of the heating/cooling wall surfaces.
     . Specify either Q (heat into the device) or Outlet Temp. 

*/
  replaceable package Medium = PartialMixtureMedium ;
  FluidPort 	inlet (redeclare package Medium = Medium)  ; 
  FluidPort 	outlet (redeclare package Medium = Medium)  ; 

  parameter Real	cf	= 1.0     ; // Pressure Loss Coeff. (m^3/Pa^0.5))     
 // parameter Area	A_we	= 1.0     ; // Wall to external external heat transfer area
  parameter Area	A_wf	= 1.0	  ; // Wall to fluid heat transfer area
//  parameter CoefficientOfHeatTransfer	h_ew	= 15	; // External to Wall heat transfer coeff. (W/m2K) 
  parameter CoefficientOfHeatTransfer	h_wf	= 150	; // Wall to fluid heat transfer coeff. 
  parameter Mass		        w_m	= 1.0	; // Mass of heat transfer walls
  parameter SpecificHeatCapacity        w_cp    = 420	; // Specific heat of wall material
  parameter Volume			holdup  = 1e-3  ; // Heater fluid holdup 
  Medium.ThermodynamicState	state 		    ; // Fluid state in the heater
  Heat			        Q_ew	    	    ; // Heat input to the device
  Heat				Q_wf		    ; // Heat tranfer from wall to fluid
  Temperature			Tw		    ; // Wall temperature (K)
  Temperature			Tf		    ; // Heater Fluid temperature (K)
  Energy			U		    ; // Internal energy of fluid holdup
  SpecificHeatCapacity          Cp    		    ; // Specific heat of fluid in holdup

  equation
    // Mass balance 
     inlet.m_flow + outlet.m_flow = 0  ;

    // Ignoring Composition change dynamics due to hold up
     outlet.Xi_outflow = inStream(inlet.Xi_outflow) ;  // Normal flow
     inlet.Xi_outflow = inStream(outlet.Xi_outflow) ;  // for  Reverse flow
     inlet.h_outflow = outlet.h_outflow ;   // Well mixed condition 

    /* Pressure differential drives flow 
       See valve eqn. for explanation  */
       inlet.m_flow =  cf * sqrt(Medium.density(state) * inlet.p)
                          * homotopy( simplified = (1 - outlet.p / inlet.p), 
                                      actual     = regRoot(1 - outlet.p / inlet.p) 
                                    );
       state = Medium.setState_pTX((inlet.p + outlet.p)*0.5 , Tf, outlet.Xi_outflow);
       outlet.h_outflow = Medium.specificEnthalpy(state);

    // Energy Balance
       inlet.m_flow * actualStream(inlet.h_outflow) 
             + outlet.m_flow *  actualStream(outlet.h_outflow) + Q_wf 
        =  der(U)  ;   // fluid
        U = holdup * Medium.density(state) * Medium.specificInternalEnergy(state) ;
	Q_wf = h_wf * A_wf * (Tw - Tf)  ;
        w_m * w_cp * der(Tw) = Q_ew - Q_wf ;  //  Wall
        Cp = Medium.specificHeatCapacityCp(state) ;
end HeaterCooler;
