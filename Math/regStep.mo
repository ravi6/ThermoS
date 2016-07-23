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

protected
  Real Alpha ;
  constant Integer SF = 5  ;

algorithm
  Alpha := 1.0 / ( 1.0 + exp(-SF * x/Tol) ) ;
  y := (1 - Alpha) * yminus + Alpha * yplus ;
//annotation(Inline="true", smoothorder=3) ;
end regStep;
