within ThermoS.Uops;
model Pipe

  // A simple  Pipe
  replaceable package Medium = PartialMixtureMedium;

  // Specify that our Medium is used in/outlet
  FluidPort inlet(redeclare package Medium = Medium);
  FluidPort outlet(redeclare package Medium = Medium);


  Medium.ThermodynamicState  state, state1, state2 ;

  parameter Length  len = 1.0  ;
  parameter Length  dia = 0.01 ;

  Real  Re, csArea, vel, tau, ff  ;

equation


  outlet.h_outflow  = inStream(inlet.h_outflow); // flow is  inlet to outlet
  inlet.h_outflow   = inStream(outlet.h_outflow); // reverse flow case

  inlet.Xi_outflow  = inStream(outlet.Xi_outflow); // pass composition un_altered
  outlet.Xi_outflow = inStream(inlet.Xi_outflow); // pass composition un_altered

  // Mass balance
  inlet.m_flow + outlet.m_flow = 0;

  // Setting the state
     state1 = Medium.setState_phX(inlet.p,  inStream(inlet.h_outflow), 
                                   inStream(inlet.Xi_outflow));
     state2 = Medium.setState_phX(outlet.p, inStream(outlet.h_outflow), 
                                   inStream(outlet.Xi_outflow));
     state = Medium.setSmoothState( (inlet.p - outlet.p), state1, state2, 100);  


     csArea =  3.147 * dia * dia / 4.0 ;
     vel   = abs(inlet.m_flow) / (csArea * Medium.density(state));  
     tau   = csArea * len / vel ;
     Re    = dia * vel *  Medium.density(state) / Medium.dynamicViscosity(state);
     ff    = if (Re < 2000) then
                  16.0 / (Re + 1e-3)   // Laminar - removing singularity at Re=0
             else
                  0.316 *  Re^(-0.25) ;  // Turbulent - Blasius Correlation
     
     inlet.p - outlet.p = sign(inlet.m_flow) *
             2 * ff * len * vel * vel  * Medium.density(state) / dia ;
end Pipe;
