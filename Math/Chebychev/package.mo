within ThermoS.Math ;
package Chebychev    "Package to compute chebychev polynomial values, derivatives and roots"  

/*
Author:  Ravi Saripalli
Version: 1.0
Date:    6th March 2015
*/


function cheby

/* Retruns Chebychev polynomial value of degree n
    at a given x    -1<x<1
*/

    input Integer n = 0;
    input Real x = 0;
    output Real y;

    algorithm
           y := cos(n * acos(x));
end cheby;

function Tnots
/* Retrun roots of a Chebychev polynomial 
     degree n 
*/
  import Modelica.Constants.pi;

    input Integer n = 0;
    output Real[n] y;

    algorithm
           y := zeros(n);  // initialize array
           for k in 1:n loop
              y[k] := cos((2 * k - 1) * pi / ( 2 * n));
           end for;
end Tnots;



function T
/* Retruns a vector of Chebychev polynomial value from
     degree 0 to n-1 a given x    -1<x<1
*/

    input Integer n = 0;
    input Real x = 0 ;
    output Real [n] y ;

    algorithm
           y := zeros(n);  // initialize array
           y[1] := 1;   y[2] :=  x;

           for k in 2:n-1 loop
              y[k+1] := 2 * x * y[k] - y[k-1];
           end for;
end T;


function Tx
/* Retruns a vector of Chebychev polynomial  
     first derivative value from
     degree 0 to n-1 a given x    -1<x<1
*/

    input Integer n = 0;
    input Real x = 0;
    output Real [n] yder;

    protected Real [n] y;

    algorithm
           y := T(n,x);  //  Chebchev values
           yder := zeros(n);  // first derivatives
           yder[1] := 0;   yder[2] :=  1;

           for k in 2:n-1 loop
              yder[k+1] := 2 * x * yder[k] - yder[k-1]  + 2 * y[k];
           end for;
end Tx;


function  Txx
/* Retruns a vector of Chebychev polynomial  
      second derivative values from
     degree 0 to n-1 a given x    -1<x<1
*/

    input Integer  n = 0;
    input Real  x = 0;
    output Real[n] yder2;

    protected Real [n] yder;

    algorithm
           yder  := Tx(n,x);  // Get array of Chebychev
           yder2 := zeros(n);  // initialize array
           yder2[1] := 0;   yder2[2] :=  0;

           for k in 2:n-1 loop
              yder2[k+1] := 2 * x * yder2[k] - yder2[k-1]  + 4 * yder[k];
           end for;

end Txx;

// ==========================================
// Functions for Shifted Chebychev Polynomial
// ==========================================

function sT
/* Retruns a vector of shifted Chebychev polynomial value from
     degree 0 to n-1 a given x    0<x<1
*/

    input Integer n = 0;
    input Real x = 0 ;
    output Real [n] y ;

    algorithm
           y := T(n, (2 * x - 1)) ;     // Reuse of unshifted function
end sT;


function sTx
/* Retruns a vector of shifted Chebychev polynomial  
     first derivative value from
     degree 0 to n-1 a given x    0<x<1
*/

    input Integer n = 0;
    input Real x = 0;
    output Real [n] yder;

    algorithm
         yder := 2 * Tx(n,  (2 * x - 1)) ;  // reuse unshifted function

end sTx;


function  sTxx
/* Retruns a vector of Chebychev polynomial  
      second derivative values from
     degree 0 to n-1 a given x    0<x<1
*/

    input Integer  n = 0;
    input Real  x = 0;
    output Real[n] yder2;

    algorithm
         yder2 := 4 * Txx(n,  (2 * x - 1)) ;  // reuse unshifted function

end sTxx;

function sTnots
/* Retrun roots of a shifted Chebychev polynomial 
     degree n 
*/
  import Modelica.Constants.pi;

    input Integer n = 0;
    output Real[n] y;

    algorithm
        y := ( Tnots(n) .+ 1 ) / 2 ;
        
end sTnots;

function sTi
/* Returns vector of integrals [0 to 1] of shifted Chebychev polynomials
   from order 0 to n-1  
*/

    input Integer  n = 0;
    output Real[n] y ;   // integrals

algorithm
        y[1] := 1.0  ; y[2] := 0.0;
        for k in 2:n-1 loop
          if (mod(k,2) == 0) then
             y[k+1] := 1.0 / (1 - k * k) ;
          else
             y[k+1] := 0 ;
          end if ;
          //   print("i="+ String(k) + " : " + String(y[k+1]) + "\n");
        end for;
      
end sTi;


//  Some miscellaneous stuff here

function cheby_brute
/* Retruns a vector of Chebychev polynomial value from
     degree 0 to n-1 a given x    -1<x<1
     Algorithm is not optimized for speed
*/

input Integer  n = 0;
input Real  x = 0;
output Real[n] y;


algorithm
       y := zeros(n);  // initialize array
       for k in 0:n-1 loop
          y[k+1] := cos(k * acos(x));
         print("i="+ String(k) + " : " + String(y[k+1]) + "\n");
       end for;
end cheby_brute;

end Chebychev  ;
