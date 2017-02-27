within ThermoS.Math;
// Date: 3rd Nov. 2016
// Author: Ravi Saripalli

function sMin
/*  A smooth Min(x,y) function, C1 continuous
    Courtesy of GAMS
*/

input Real x        ;
input Real y        ;
input Real eps      ;  // Tolerance of error
output Real z       ;


algorithm
 z := - ( sqrt((x - y)^2 + eps^2 ) - (x - y) ) / 2 ;
annotation(Inline="true", smoothOrder=2) ;
end sMin ; 
