within ThermoS.Uops.Tanks;

model SimpleTank    " An Open Tank "
   

/*It has two connections
     inlet ... 
     outlet ... 
*/

  replaceable   package Medium = PartialMixtureMedium ;
  FluidPort 	inlet (redeclare package Medium = Medium)  ; 
  FluidPort 	outlet (redeclare package Medium = Medium)  ; 


  Heat				Q_in		    ; // Heat tranferred to fluid  
  Temperature			Tf		    ; // Tank Fluid temperature (K)
  Pressure                      Pa                  ; // Ambient Pressure
  Energy			U		    ; // Internal energy of fluid fvol
  SpecificHeatCapacity          Cp    		    ; // Specific heat of fluid in fvol
  Medium.ThermodynamicState	state 		    ; // Fluid state in the tank
  Mass                          mf                  ; // Mass of fluid in the tank
  Volume                        fvol                ; // Fluid Volume
  Volume                        tvol                ; // Tank Volume

// Tank State
  Density                       rho                 ; // Fluid Density

  equation
    state = Medium.setState_pTX(Pa , Tf, outlet.Xi_outflow);
    outlet.h_outflow = Medium.specificEnthalpy(state);
    rho = Medium.density(state) ;

    der(mf) =  smooth(1, if(fvol < tvol) then
                            (inlet.m_flow + outlet.m_flow)    
                         else
                            0  );    // Over flow condition
    fvol = mf / rho ;

    // Energy Balance
       inlet.m_flow * actualStream(inlet.h_outflow) 
             + outlet.m_flow *  actualStream(outlet.h_outflow) + Q_in 
        =  der(U)  ;   // fluid
        U = fvol * Medium.density(state) * Medium.specificInternalEnergy(state) ;
        Cp = Medium.specificHeatCapacityCp(state) ;

    // Ignoring Composition change dynamics due to hold up
     outlet.Xi_outflow = inStream(inlet.Xi_outflow) ;  // Normal flow
     inlet.Xi_outflow = inStream(outlet.Xi_outflow) ;  // for  Reverse flow
     inlet.h_outflow = outlet.h_outflow ;   // Well mixed condition 


end SimpleTank;
