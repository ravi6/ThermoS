within ThermoS.Uops;
model OnOff
 /* An On/Off Controller with dead band */
  parameter Real mvMin = 0;
  parameter Real mvMax = 100;
  parameter Real pvMin = 0;
  parameter Real pvMax = 100;
  parameter Real deadBand = 5 ;
                             
  Real sp(min = pvMin, max = pvMax); // Setpoint
  Real pv; //  Present value 
  Real mv(min = mvMin, max = mvMax); // Controller Output (manipulated)
  Real err; // Setpoint - Present value 

equation
  err = (sp - pv) ;
  // On/Off with dead band
  mv = noEvent(if pv >= (sp + deadBand) then 
                  mvMin 
               elseif pv <= (sp - deadBand)  then 
                  mvMax
               elseif (der(pv) > 0)  then
                  mvMax
               elseif (der(pv) < 0)  then
                  mvMin
               else
                  mvMin
               );
initial equation
end OnOff;
