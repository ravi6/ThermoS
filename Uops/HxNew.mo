// Version 1.3
//  24th Mar. 2014
within ThermoS.Uops;
model HxNew

  /* This heater/cooler Model 
     . Fluid is completely mixed in the heater (CSTR)
     . Permits pressure loss specification through loss coeff
       based on volumetric flow, or as a fixed value.
     . Accounts for  thermal inertia of the heating/cooling wall surfaces.
     . Specify either Q (heat into the device) or Outlet Temp. 
     . No reverse flow

*/
  replaceable package Medium = PartialMixtureMedium ;
  FluidPort  inlet (redeclare package Medium = Medium, m_flow (max = 0))  ; 
  FluidPort outlet (redeclare package Medium = Medium, m_flow (min = 0)) ; 

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
  SpecificHeatCapacity          Cp    		    ; // Specific heat of fluid in holdup
  Volume                        holdup              ;

  equation

    // Holdup Conditions
       medium.Xi = inlet.Xi_outflow ;
       medium.h = outlet.h_outflow ; 
       medium.p = (inlet.p + outlet.p) / 2 ; 
       Tf = medium.T  ;

    // No mass accumulation
       inlet.m_flow + outlet.m_flow = 0 ;

    // No change in composition
       inlet.Xi_outflow = outlet.Xi_outflow ;
       outlet.Xi_outflow = inStream(inlet.Xi_outflow) ;

    /* Pressure differential drives flow 
       See valve eqn. for explanation  */
       inlet.m_flow =  cf * sqrt(medium.d * inlet.p)
                          * homotopy( simplified = (1 - outlet.p / inlet.p), 
                                      actual = regRoot (1 - outlet.p / inlet.p) 
                                    );
    // Energy Balance (ignore vdp/dt) 
        inlet.h_outflow = inStream (inlet.h_outflow) ;
        holdup * der ( medium.d * outlet.h_outflow) = 
                  inlet.m_flow * (inlet.h_outflow) 
                    + outlet.m_flow *  outlet.h_outflow + Q_wf  ;

	Q_wf = h_wf * A_wf * (Tw - Tf)  ;
        w_m * w_cp * der(Tw) = Q_ew - Q_wf ;  //  Wall
        Cp = Medium.specificHeatCapacityCp(medium.state) ;

        //assert(inlet.m_flow <= 1e-6, "Reverse flow at inlet!");
        //assert(outlet.m_flow >= 1e-6, "Reverse flow at outlet!");

end HxNew;
