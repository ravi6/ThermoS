within ThermoS.Uops;
model Controller
 /* A PID Controller  Controller */
  Real 	Kc  ; // Proportional Gain 
  Real 	Ti  ; // Integral Time
  Real	Td  ; // Derivative Time
  parameter Boolean 	reverseActing = true;
                        /* Manipulated varibale decreases when
                             measured variable increases */
  parameter Real mvMin = 0;
  parameter Real mvMax = 100;
  parameter Real pvMin = 0;
  parameter Real pvMax = 100;
  parameter Integer action = if reverseActing then -1
                             else 1 ;
                             
  Real sp(min = pvMin, max = pvMax); // Setpoint
  Real pv; //  Present value 
  Real mv(min = mvMin, max = mvMax); // Controller Output (manipulated)
  Real err; // Setpoint - Present value 
  Real intErr; // error integral
equation
  err = (sp - pv) ;
  // Integral Action with (No windup)
  der(intErr) = noEvent(if mv < mvMin and err < 0 
                                  or 
                           mv > mvMax and err > 0 
                        then 0 else err);
  mv  = action * ((mvMax - mvMin) / (pvMax - pvMin)) 
         *   Kc * ( err  + intErr / Ti  + Td * der(err) )  + mvMin ;   
initial equation
  intErr = 0;
end Controller;
