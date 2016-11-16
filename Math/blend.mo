within ThermoS.Math;
// Date: 10th Nov 2016
// Author: Ravi Saripalli

/* Facilitate blending of step changes
 When  x >> Tol  y -> 1 
         x = 0     y = 0.5
         x << -Tol  y -> 0
  x is assumed to be appropriately scaled
*/

model blend
 parameter Real Tol=0.02              ;
 Real x (min=-20, max=20, start=0)    ;
 Real f (min=0, max=1, start=0.5)     ;

equation
  f =   1.0 / (1.0 + exp(-x/Tol))  ;
end blend;
