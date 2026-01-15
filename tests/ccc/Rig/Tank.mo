within Rig;
model Tank
/*
*  A Gas Storage Vessel 
*/

    replaceable package Medium = PartialMixtureMedium ;
    FluidPort inlet(redeclare package Medium = Medium) ;
    FluidPort outlet(redeclare package Medium = Medium) ;  

//  Parameters
  parameter    Volume   	vol = 1    ;   // Tank Volume (m3)
  parameter    Real	        cf  = 1.0  ; // Pressure Loss Coeff. (m^3/Pa^0.5)) 

// State Variables
    Medium.Temperature		T		;
    Medium.AbsolutePressure	p		;
    EnthalpyFlowRate            Q_loss          ; // Heatloss from tank

    Medium.BaseProperties       medium ;
    Mass                 M  ; // Gas mass in tank
    InternalEnergy       U  ; // of tank contents
    
  equation

    // Tank p,T for ease of access
    T = medium.T  ;
    p = medium.p  ;
    medium.Xi = Medium.reference_X [1:Medium.nXi] ;
    
    // Mass Balance
    M = medium.d * vol ;
    der(M) = inlet.m_flow + outlet.m_flow ;

     // Enthalpy Balance
     U = M * medium.u ;
     der(U)  = - Q_loss 
                     + inlet.m_flow * actualStream(inlet.h_outflow)
                     + outlet.m_flow * actualStream(outlet.h_outflow);

     // Well Mixed Condition
      outlet.Xi_outflow = inlet.Xi_outflow ;
      inlet.Xi_outflow = medium.Xi ;
      inlet.h_outflow  = outlet.h_outflow  ;
      outlet.h_outflow = medium.h ;

     // Inlet and outlet port resistences (assume to be same)
       inlet.m_flow =  cf * sqrt(medium.d * inlet.p)
                          * homotopy( simplified = (1 - medium.p / inlet.p), 
                                      actual     = regRoot(1 - medium.p / inlet.p) 
                                    );

       outlet.m_flow =  - cf * sqrt(medium.d * medium.p)
                          * homotopy( simplified = (1 - outlet.p / medium.p), 
                                      actual     = regRoot(1 - outlet.p / medium.p) 
                                    );
end Tank;
