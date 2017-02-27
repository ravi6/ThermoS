within ThermoS.Math;
model filter

// A firstorder filter 

Real y ;
Real x ;
parameter Real tau = 1 ;

equation
 tau *  der(y) +  y = x ;

initial equation
 y = x ; 

end filter;
