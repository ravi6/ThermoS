within Rig;

model HeatX

  /* This heater/cooler Model 
     . Fluid is completely mixed in the heater (CSTR)
     . Accounts for  thermal inertia of the heating/cooling wall surfaces.
     . Specify either Q (heat into the device) or Outlet Temp. 
     . No reverse flow permitted.
*/
  replaceable package Medium = PartialMixtureMedium ;
  FluidPort  inlet (redeclare package Medium = Medium, m_flow (min = 0))  ; 
  FluidPort outlet (redeclare package Medium = Medium, m_flow (max = 0)) ; 

  parameter Real	cf	= 1.0     ; // Pressure Loss Coeff. (m^3/Pa^0.5)) 
  parameter Volume      holdup = 0.1      ;
  parameter Area	A_wf	= 1.0	  ; // Wall to fluid heat transfer area
  parameter CoefficientOfHeatTransfer	h_wf	= 150	; // Wall to fluid heat transfer coeff. 
  parameter Mass		        w_m	= 1.0	; // Mass of heat transfer walls
  parameter SpecificHeatCapacity        w_cp    = 420	; // Specific heat of wall material

  Medium.BaseProperties         medium  ;
  HeatFlowRate		        Q_ew	; // Heat input to the device
  HeatFlowRate			Q_wf	; // Heat tranfer from wall to fluid
  Medium.Temperature		Tw	; // Wall temperature (K)
  Medium.Temperature		Tf	; // Heater Outlet Fluid temperature (K)
  Medium.AbsolutePressure	p	; // Heater fluid pressure (Pa)


  Mass                 M       ; // Mass of holdup
  Enthalpy             H       ; // Enthalpy of holdup

  equation

     // Mass Balance
       M = holdup * medium.d  ;
       der (M) = inlet.m_flow + outlet.m_flow  ;  


    // Energy Balance
	Q_wf = h_wf * A_wf * (Tw - Tf)  ;  // Wall to Fluid Heat transfer 
        w_m * w_cp * der(Tw) = Q_ew - Q_wf ;  //  Thermal Inertia of Wall
    
        H = M * medium.h ;
        der (H) =  inlet.m_flow * actualStream (inlet.h_outflow)
                   + outlet.m_flow * actualStream (outlet.h_outflow)
                   + Q_wf  + holdup * der (p) ;

     // State of Holdup
        medium.Xi = Medium.reference_X [1:Medium.nXi] ;
        medium.h = outlet.h_outflow ;
        medium.T = Tf ;   
        medium.p = p  ;

     // Well Mixed Condition
        inlet.h_outflow = outlet.h_outflow ; 
        inlet.Xi_outflow = outlet.Xi_outflow ;
        outlet.Xi_outflow = medium.Xi ;

     // Inlet and outlet port resistences (assume to be same)
       inlet.m_flow =  cf * sqrt(medium.d * inlet.p)
                          * homotopy( simplified = (1 - medium.p / inlet.p), 
                                      actual     = regRoot(1 - medium.p / inlet.p) 
                                    );

       outlet.m_flow =  - cf * sqrt(medium.d * medium.p)
                          * homotopy( simplified = (1 - outlet.p / medium.p), 
                                      actual     = regRoot(1 - outlet.p / medium.p) 
                                    );
//        assert(inlet.m_flow >= 1e-6, "Reverse flow at inlet!");
//        assert(outlet.m_flow <= 1e-6, "Reverse flow at outlet!");

      
end HeatX;
