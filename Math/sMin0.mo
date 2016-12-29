within ThermoS.Math;
// Date: 3rd Nov. 2016
// Author: Ravi Saripalli

function sMin0
/*  A smooth Min(x,0) function, C1 continuous
    Courtesy of GAMS
*/

input Real x        ;
input Real eps      ;  // Tolerance of error
output Real y       ;


algorithm
 y := - ( sqrt( (x * x) + (eps * eps) )  - x  ) / 2 ;
annotation(Inline="true", smoothOrder=2) ;
end sMin0 ; 
