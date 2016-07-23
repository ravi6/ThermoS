within ThermoS.Uops.Valves;

/* This Valve has some time constant (first order)
      It can't change its position instantaneously
*/
model RealValve                        // A control Valve
    extends ThermoS.Uops.Valves.Valve  ;

parameter Time tau = 0.01 ; // 10 milliseconds (default)
Percent   spo   ; // Desired valve openning 

equation
      tau * der(po) + po = spo ; 

end RealValve;
