// Version 1.3
//  24th Mar. 2014
within ThermoS.Uops.HeatExch;
model HeatX

  /* This heater/cooler Model 
     . Fluid is completely mixed in the heater (CSTR)
     . Permits pressure loss specification through loss coeff
       based on volumetric flow, or as a fixed value.
     . Accounts for  thermal inertia of the heating/cooling wall surfaces.
     . Specify either Q (heat into the device) or Outlet Temp. 
     . No reverse flow permitted.
     . Symbolic derivative friendly version
*/
  replaceable package Medium = PartialMixtureMedium ;
  FluidPort  inlet (redeclare package Medium = Medium, m_flow (min = 0))  ; 
  FluidPort outlet (redeclare package Medium = Medium, m_flow (max = 0)) ; 

  parameter Real	cf	= 1.0     ; // Pressure Loss Coeff. (m^3/Pa^0.5)) 
  parameter Area	A_wf	= 1.0	  ; // Wall to fluid heat transfer area
  parameter CoefficientOfHeatTransfer	h_wf	= 150	; // Wall to fluid heat transfer coeff. 
  parameter Mass		        w_m	= 1.0	; // Mass of heat transfer walls
  parameter SpecificHeatCapacity        w_cp    = 420	; // Specific heat of wall material

  Medium.BaseProperties         medium              ;
  HeatFlowRate		        Q_ew	    	    ; // Heat input to the device
  HeatFlowRate			Q_wf		    ; // Heat tranfer from wall to fluid
  Medium.Temperature		Tw		    ; // Wall temperature (K)
  Medium.Temperature		Tf		    ; // Heater Fluid temperature (K)
  Volume                        holdup ;

  equation

    // Holdup Conditions
       medium.p = (inlet.p + outlet.p) / 2 ;
       medium.Xi = outlet.Xi_outflow ;
       medium.T = Tf ;   // For ease of access to outlet temp

    // Mass balance 
       holdup * der (medium.d) = inlet.m_flow + outlet.m_flow;   

    // Ignoring Composition change dynamics due to hold up
    // Holdup is well Mixed
       outlet.Xi_outflow = inStream(inlet.Xi_outflow) ;  // Normal flow
       inlet.Xi_outflow = outlet.Xi_outflow ;
       inlet.h_outflow = outlet.h_outflow ;

    /* Pressure differential drives flow 
       See valve eqn. for explanation  */
       inlet.m_flow =  cf * sqrt(medium.d * inlet.p)
                          * homotopy( simplified = (1 - outlet.p / inlet.p), 
                                      actual     = regRoot(1 - outlet.p / inlet.p) 
                                    );
    // Energy Balance
       holdup * der(medium.d * medium.h)   
           = inlet.m_flow * actualStream (inlet.h_outflow) 
             + outlet.m_flow * actualStream (outlet.h_outflow) 
             + holdup * der (medium.p) + Q_wf ;

	Q_wf = h_wf * A_wf * (Tw - Tf)  ;
        w_m * w_cp * der(Tw) = Q_ew - Q_wf ;  //  Wall

//        assert(inlet.m_flow >= 1e-6, "Reverse flow at inlet!");
//        assert(outlet.m_flow <= 1e-6, "Reverse flow at outlet!");

end HeatX;
