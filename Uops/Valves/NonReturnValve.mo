within ThermoS.Uops.Valves;

model NonReturnValve
    extends ThermoS.Uops.Valves.partialValve  ;

equation

    inlet.m_flow = smooth(0,
                          if (inlet.p > outlet.p) then
                             cv * sqrt(Medium.density(state)) * sqrt(abs(inlet.p - outlet.p))
                          else
                             0 
                          );

end NonReturnValve;
