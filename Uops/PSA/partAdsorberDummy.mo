within ThermoS.Uops.PSA ;

partial model partAdsorberDummy
import ThermoS.Math.Chebychev.*  ;

/*  Adsorption Column Model 
                  Author: R. Saripalli
                  Date:   10th March 2015 
                  Last Modified:  12th Jan 2016
*/

/*  Extend this model to provide necessary B.Cs and I.Cs
       bedParams Record contains all the necessary input data
       note: the default Record contains bed data for Zeolite X5
             for N2/O2 separation
*/

parameter Real uTol = 1e-5 ;  // b.c change regularized Step Tolerence

constant Integer Nc = 2 ;
constant Integer N =  10 ;

// Chebychev package functions used to get all the data needed for collocation
constant Real [:] zi = sTnots(N-2)        ;  // Interior Collocation Points
constant Real [:] z  = cat(1, zi, {0, 1}) ;  // Including boundaries

constant Real [:,:] vT   = transpose(sT(N, z))     ;  // Cheby values (index i) at z (index j)
constant Real [:,:] vTz  = transpose(sTx(N, z))    ;  // first deriv Cheby values (index i) at z (index j)
constant Real [:,:] vTzz = transpose(sTxx(N, z))   ;  // second deriv Cheby values (index i) at z (index j)
constant Real [:]   vTi  = sTi(N)                  ;  // Cheby integrals over 0 to 1

// Add some limits to variables (might assist in convergence ???)
type Frac = Real(min=0, max=1, start=0.1, nominal=1);
type Vel  = Real(min=-130, max=130, nominal = 0.1, start = 0.01) ;
type Conc = Real(min=0, max=100, nominal = 1, start = 1.0); 
type Coef = Real(min=-120, max=120) ;  // Don't bad guess these lest you have problems
type dCoef = Real(min=-120, max=120) ;  // time Derivatives of Coefs
type Src  = Real(min=-200, max=200); // , start = 0); 
type Press = Real(min=1e-2, max=20, start=1, nominal=1);
type Temp = Real(min=0.9, max=5, start=1, nominal=1); 

Coef [N, Nc-1] Coef_y   ;    // Coeffs. of collocation for gas mole fraction
Coef [N]       Coef_p   ;    // Coeffs. of collocation for gas pressure
Coef [N]       Coef_u   ;    // Coeffs. of collocation for gas pressure
Coef [N, Nc]   Coef_Q   ;    // Coeffs. of collocation for adsorbed concentration

// Time derivatives of Coefs. (These become dummyDerivatives ... Reducing the index)
//    in otherwords, we will transfrom our ODE's set based on Coefs to actual variables
//    these coefs. relate to. Naturally we need to add that many more eqns. But it is
//    worth it I think.

dCoef [N, Nc-1] dCoef_y   ;    // time Derivatives of Coeffs. of collocation for gas mole fraction
dCoef [N]       dCoef_p   ;    // time Derivatives of Coeffs. of collocation for gas pressure
dCoef [N, Nc]   dCoef_Q   ;    // time Derivatives of Coeffs. of collocation for adsorbed concentration

Frac   [N, Nc]   y      ;    // Gas mole fractions at collocation points
Vel    [N]       u      ;   // Gas velocity at collocation pionts
Press  [N]       p      ;   // Gas pressure at collocation pionts


Conc [N, Nc]   Q        ;   // Adsorbed concentraion in solid  at collocation points
Conc [N, Nc]   Qeq      ;   // Adsorbed eqilibrium concentraion in solid  at collocation points
Src  [N, Nc]   S        ;   // Adsorbed rate  at collocation points

// Boundary Variables
Frac   [Nc]      yin_in, yin_out       ;    // Gas mole fractions (z=0- & z=1+ used when flow is into bed)
Press            p_in, p_out       ;    // Gas pressure at upstream (z=0, z=1)


Real [N] zs   ;  // Including boundaries

equation

zs = z ;       // I need this as I can't access constants in post processing


// Nodal Values of pressure and velocity
p =    Coef_p * vT ;
u =    - (1.0 / bedParams.Kappa) * (Coef_p * vTz) ;  // pressure gradient determines velocity
u = Coef_u * vT ;  // makes my plotting easy  (this is just curve fitting)


// DummyDerivatives (Time derivatives of Coefs. in our case)
// Now we see time derivatives of real variables, and avoid Coefs being states

der(p) = dCoef_p * vT ;

for n in 1:Nc-1 loop              
   der(y[:, n]) =   dCoef_y[:, n] * vT  ;  
end for ;

for n in 1:Nc loop              
   der(Q[:, n]) =   dCoef_Q[:, n] * vT  ; 
end for ;


// Nodal Values of Gas phase compostion, adsorbed/equil Concentrations in solid
// and adsoption rate into solid

for n in 1:Nc-1 loop              
      y[:, n] =   Coef_y[:, n] * vT    ;
end for ;

for i in 1:N loop
       sum ( y[i, n] for n in 1 : Nc ) = 1 ; 
end for ; 

for n in 1:Nc loop              
      Q[:, n] =  Coef_Q[:, n] * vT   ;
    for m in 1:N loop
           Qeq[m, n] = max(0, 
                           bedParams.Qs[n] * ( bedParams.B[n] *  p[m] * y[m, n] ) 
                           / ( bedParams.Tb + p[m] * sum ( bedParams.B[j] *  y[m, j] for j in 1 : Nc) )
                          ) ;
          S[m, n]    =    bedParams.Km[n] * ( Qeq[m, n] - max(0, Q[m, n]) );  // Rate of adsorption
    end for ;
end for ;


// Component Balance Eqns.
 for  n in  1:Nc-1 loop       // nth Component Balance
     for m in 1:N-2 loop      // interior Collocation
         p[m] * ( dCoef_y[:, n] * vT[:, m]                               // p * (dy_i/dt
                 +  u[m]  *  Coef_y[:, n] * vTz[:, m] )                  //       u * dy_i/dz)
                 - (1.0 / bedParams.Pe) * (                              //  - (1/Pe)*( 
                    p[m] * Coef_y[:, n] * vTzz[:, m]                     //       p * d2y_i/dz2 
                    + y[m, n] * (  Coef_p * vTzz[:, m]                   //       + yi * (d2p/dz2 
                                    + 2 * (Coef_p * vTz[:, m])           //                 + 2 * dp/dz 
                                        * (Coef_y[:, n] * vTz[:, m]) ))  //                  * dyi/dz))
             +  ( S[m, n]                                                // ( dQ_i/dt
                  - y[m, n]                                              //    - y_i
                    * sum ( S[m, j] for j in 1 : Nc )                    //         * sum (dQ_i/dt)
                ) * bedParams.Epsilon * bedParams.Tb = 0 ;               //   ) * epsilon * T 

     end for; //end of collocation point
end for; // end of all component balances

// Total Mass Balance Equaiton 
     for m in 1:N-2 loop      // interior Collocation
       dCoef_p * vT[:, m]                      // dp/dt  
            + u[m] * Coef_p * vTz[:, m]        // + u * dp/dz
            + p[m] * ( - (1.0 / bedParams.Kappa) * Coef_p * vTzz[:, m] )        // + p * du/dz
               + bedParams.Epsilon                                      // (1-e)/e
                    *  bedParams.Tb * sum ( S[m, j] for j in 1 : Nc )   //    T     * sum (dQ_i/dt)
               = 0 ; 
     end for;


// Adsorption rate into the solid
 for  n in  1:Nc loop       // nth Component 
   for m in 1:N loop      // interior Collocation
        dCoef_Q[:, n] * vT[:, m]                             // dQ_i/dt
                         - S[m, n]                           // - Km_i * (Qi_e - Qi)
          = 0 ;   
     end for;
 end for; // end of all component balances

// Boundary Conditions

 for  n in  1:Nc-1 loop       // nth Component 
    Coef_y[:, n] * vTz[:, N-1] =  regStep(u[N-1], // (p_in - p[1]), 
                                      bedParams.Pe * u[N-1] * (y[N-1, n] - yin_in[n]), 0,
                                      uTol)   ; // bed inlet

    Coef_y[:, n] * vTz[:, N]   =  regStep(u[N], // (p[N] - p_out), 
                                    0, bedParams.Pe * u[N] * (y[N, n] - yin_out[n]),
                                    uTol)  ; // bed outlet 
 end for;

// Pressure Boundary Conditions
     Coef_p[:] * vT[:, N-1] = p_in ;    // Inlet  pressure 
     Coef_p[:] * vT[:, N] = p_out ;    // outlet  pressure 
//     Coef_p[:] * vTz[:, N]   = 0 ;   // outlet pressure  (zero gradient = no flow)

end partAdsorberDummy;


      // y[i, Nc] = 1 .- min(1, sum ( y[i, n] for n in 1 : Nc-1 )) ; 
      //y[i, Nc] = 1 .-  sum ( y[i, n] for n in 1 : Nc-1 ) ; 
