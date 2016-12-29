model fo
Real y ;
Real x ;
parameter Real tau = 10 ;

equation
 tau *  der(y) +  y = x ;

initial equation
 y = x ; 

end fo;
