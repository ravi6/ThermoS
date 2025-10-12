within ThermoS.Uops;
model CompressorBasic
  "Polytropic compression model (algebraic, steady flow)"

  replaceable package Medium = PartialMixtureMedium;

  // Ports
  Modelica.Fluid.Interfaces.FluidPort_a inlet(redeclare package Medium = Medium);
  Modelica.Fluid.Interfaces.FluidPort_b outlet(redeclare package Medium = Medium);

  // Parameters
  parameter Real n(min = 1, max = 1.5) = 1.4 "Polytropic exponent";
  parameter Real pr(min = 0.1, max = 10, start = 2) "Pressure ratio (p_out/p_in)";
  parameter Real effPoly(min = 0, max = 1) = 0.9 "Polytropic efficiency";
  parameter Boolean adiabatic = true "If true, Q = 0";
  parameter Power Ws = 1000 "Shaft power input (positive to fluid)";
  parameter HeatFlowRate Q = 0 "Heat transfer to fluid (+ve in)";

  // Medium states
  Medium.ThermodynamicState state_in;
  Medium.ThermodynamicState state_out;

  // Flow properties
  Medium.MassFlowRate m_flow(start = 0.1, fixed = false);
  Medium.SpecificEnthalpy h_outflow(start = 3e5, fixed = false);

  Real gamma(start = 1.4);

equation
  // === Port definitions ===
  inlet.m_flow + outlet.m_flow = 0;
  outlet.p = pr * inlet.p;

  // === Define states ===
  state_in  = Medium.setState_phX(inlet.p,  inStream(inlet.h_outflow),
                                  inStream(inlet.Xi_outflow));
  state_out = Medium.setState_phX(outlet.p, h_outflow,
                                  inStream(inlet.Xi_outflow));

  gamma = Medium.isentropicExponent(state_in);

  // === Energy balance ===
  // Q = m*(h_out - h_in) - Ws  → relates m_flow and h_outflow
  Q = m_flow * (Medium.specificEnthalpy(state_out)
                - Medium.specificEnthalpy(state_in)) - Ws;

  // === Polytropic relation (temperature ratio) ===
  state_out.T = state_in.T * pr^((gamma - 1)/(gamma * effPoly));

  // === Port enthalpy connection ===
  outlet.h_outflow = h_outflow;
  inlet.Xi_outflow = outlet.Xi_outflow;
  outlet.Xi_outflow = inStream(inlet.Xi_outflow);

  // Adiabatic option enforcement
  if adiabatic then
    Q = 0;
  end if;

end CompressorBasic;

