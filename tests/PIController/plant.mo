model plant
 /* A PID Controller  Controller */

  parameter Real 	Kc = 2  ; // Proportional Gain (%)
  parameter Real 	Ti = 20  ; // Integral Time
  parameter Real	Td = 1e-3  ; // Derivative Time
  parameter Boolean 	reverseActing = false;
                        /* Manipulated varibale decreases when
                             measured variable increases */

  parameter Real mvMin = 0;
  parameter Real mvMax = 100;
  parameter Real pvMin = 0;
  parameter Real pvMax = 100;
  parameter Integer action = if reverseActing then -1
                             else 1 ;
                             

  Real sp; // Setpoint
  Real pv; // Measured value 
  Real mv(min = 0, max = 1); // Controller Output (manipulated value)
  Real err; // Setpoint - Measurevalue
  Real intErr; // error integral
  Real x, derx, y;
  Real BandGain ;

equation
  err = (sp - pv) ;
  // Integral Action with (No windup)
  der(intErr) = noEvent(if mv < mvMin and err < 0 
                                  or 
                           mv > mvMax and err > 0 
                        then 0 else err);

  BandGain = action * 100 * (mvMax - mvMin) / (pvMax - pvMin) ;
  mv = BandGain * Kc * ( err  + intErr / Ti  + Td * der(err) ) ;   

  derx = der(x);
  y = der(derx) + 2 * 0.5 * der(x) + x;

  sp = 0.5;
  pv = y;
  mv = x;
initial equation
  intErr = 0;
  derx = 2.0;
  
end plant;
