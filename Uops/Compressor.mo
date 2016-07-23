within ThermoS.Uops;
model Compressor

/* Version 1.5 
   Author: Ravi Saripalli
   Last Modified 13th March 2014
Change Log:
   Compressor flow is limited under choked conditions
*  A Gas Compressor with two ports
*     Note in this model
*      inlet and outlet (do not permit flow reversal)
*        Rotor dynamics, choke and stonewall conditions handled
*       Set power input to compressor speed adjusts to it.
*/

  replaceable package CompMap = TurboMachineMap;
  replaceable package Medium = PartialMixtureMedium;
  FluidPort  inlet(redeclare package Medium  = Medium, m_flow(min=0));
  FluidPort outlet(redeclare package Medium  = Medium, m_flow(max=0));
  parameter MomentOfInertia 	rotor_MI	 =1; // kg.m2

  // Diagnostic Variables
  MassFlowRate 		mflow_choke;    // MassFlowRate mflow_surge;
  Real 			beta(min=0, max=1); // Parametric Line Value at operating point
  MassFlowRate 		Mc;
  Real 			Nc;
  Real 			eff(min=0, max=1); //(adiabatic efficiency)

  // State Variables
  Medium.ThermodynamicState state_in, state_out;
  Medium.SpecificEntropy s_in, s_out;
  AngularVelocity 	speed(min=0.1); //(avoid singularity near zero speed)
  Real 			pratio(min=1);
  Power 		Ws;
  Power 		power ; // Power input to compressor

equation
  // Mass and Component Balance
  inlet.m_flow + outlet.m_flow = 0; // No accumulation

  // Just pass composition across the device
  inStream(inlet.Xi_outflow) = outlet.Xi_outflow;
  inStream(outlet.Xi_outflow) = inlet.Xi_outflow;
  inlet.h_outflow = inStream(inlet.h_outflow); // Set device inlet enthalpy

  // Energy Balance
  inlet.m_flow*actualStream(inlet.h_outflow) 
	+ outlet.m_flow*actualStream( outlet.h_outflow) + Ws = 0;

  /* Compressor Map Lookup 
        Corrected Speed and Flows */
           Mc = CompMap.fnMc(state_in.T, inlet.p, inlet.m_flow);
           Nc = CompMap.fnMc(state_in.T, inlet.p, speed);
           mflow_choke = CompMap.fnChokedFlow(speed, state_in.T, inlet.p);

  // Pressure Rise across the unit
         outlet.p = inlet.p*pratio;    
         pratio = CompMap.getPratio (beta, Nc) ;
         eff    = CompMap.getEff (beta, Nc) ;
         Mc     = CompMap.getMc (beta, Nc) ;

  // Outflow Enthalpy Calculation
  state_in = Medium.setState_phX( inlet.p, inlet.h_outflow, inlet.Xi_outflow);
  outlet.h_outflow = inlet.h_outflow +
                   (Medium.isentropicEnthalpy(outlet.p, state_in) - inlet.h_outflow)/eff;
  state_out = Medium.setState_phX( outlet.p, outlet.h_outflow, outlet.Xi_outflow);

  // Rotor Dynamics
  rotor_MI*speed*der(speed) =  power -  Ws ;

  // Debug Code
  s_in = Medium.specificEntropy(state_in);
  s_out = Medium.specificEntropy(state_out);
end Compressor;
