within ThermoS.Math;
// Date: 19th July 2016
// Author: Ravi Saripalli

function regStep
/* Blend a Step function with a Sigmoidal Shape
     when x reaches +/- Tol we get stepped out put
     and  in between we have the smooth belend wihtin +/-Tol
*/

input Real x        ;
input Real yplus    ;
input Real yminus   ;
input Real Tol      ;
output Real y       ;


algorithm
  y :=           ( 1.0 / (1.0 + exp(-x/Tol)) )  * yplus
       +  (1.0 - ( 1.0 / (1.0 + exp(-x/Tol)) )) * yminus ;
annotation(Inline="true", smoothorder=3) ;
end regStep;
