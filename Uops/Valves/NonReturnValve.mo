within ThermoS.Uops.Valves;

model NonReturnValve
    extends ThermoS.Uops.Valves.partialValve  ;

Percent      po        ;               // Valve % Open  

equation

    inlet.m_flow = noEvent(smooth(0,
                          if (inlet.p > outlet.p) then
                             cv * (po / 100.0)  * sqrt(med.d) * sqrt(abs(inlet.p - outlet.p) + 1)
                          else
                             0 
                          ));

end NonReturnValve;
