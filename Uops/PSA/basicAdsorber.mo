within ThermoS.Uops.PSA ;

model basicAdsorber  extends partAdsorber;

/*  Adsorption Column Model 
                  Author: R. Saripalli
                  Date:   10th March 2015 
                  Last Modified:  2nd Jul 2015
*/


equation

// Pressure Boundary Conditions
     p_in =  1.5 ; //+ 2 * (1 - exp(-time/0.1)) ;
     p_out = 1.0 ; //p_in;  //pressurizing from both ends

// Compostion Boundary Conditions (used when flow is into the bed)
    yin_in = {0.79, 0.21} ;
    yin_out = {0.79, 0.21} ; 

// Initial Conditions
initial equation

  for n in 1:Nc-1 loop             // Gas Composition in interior 
     for m in 1:N-2 loop
       y[m, n]  = yin_in[n]   ;
     end for ;
  end for ;

   for m in 1:N loop
    Q[m,1] = 2.9 ;    // Adsorbate concentrations
    Q[m,2] = 0.37 ;
   end for ;

   for m in 1:N-2 loop
    p[m] = 1.0 ;
   end for ;

end basicAdsorber;
