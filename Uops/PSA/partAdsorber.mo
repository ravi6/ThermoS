within ThermoS.Uops.PSA ;

partial model partAdsorber
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

parameter Real pTol = 1e-3 ;  // b.c change regularized Step Tolerence

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
type Frac = Real(min=0, max=1, start=0.5, nominal=1);
type Vel  = Real(min=-130, max=130, nominal = 0.1, start = 0) ;
type Conc = Real(min=0, max=100, nominal = 1, start = 1.0); 
type Coef = Real(min=-120, max=120) ;  // Don't bad guess these lest you have problems
type Src  = Real(min=-200, max=200, start=0); // , start = 0); 
type Press = Real(min=1e-2, max=20, start=1, nominal=1);
type Temp = Real(min=0.9, max=5, start=1, nominal=1); 

Coef [N, Nc-1] Coef_y   ;    // Coeffs. of collocation for gas mole fraction
Coef [N]       Coef_p   ;    // Coeffs. of collocation for gas pressure
Coef [N]       Coef_u   ;    // Coeffs. of collocation for gas pressure
Coef [N, Nc]   Coef_Q   ;    // Coeffs. of collocation for adsorbed concentration

Frac   [N, Nc]   y(each stateSelect=StateSelect.always)      ;    // Gas mole fractions at collocation points
Vel    [N]       u      ;   // Gas velocity at collocation pionts
Press  [N]       p(each stateSelect=StateSelect.always)     ;   // Gas pressure at collocation pionts


Conc [N, Nc]   Q(each stateSelect=StateSelect.always)        ;   // Adsorbed concentraion in solid  at collocation points
Conc [N, Nc]   Qeq      ;   // Adsorbed eqilibrium concentraion in solid  at collocation points
Src  [N, Nc]   S        ;   // Adsorbed rate  at collocation points

// Boundary Variables
Frac   [Nc]      yin_in, yin_out       ;    // Gas mole fractions (z=0- & z=1+ used when flow is into bed)
Press            p_in(start=1), p_out(start=1)       ;    // Gas pressure at upstream (z=0, z=1)


Real [N] zs   ;  // Including boundaries

equation

zs = z ;       // I need this as I can't access constants in post processing


// Nodal Values of pressure and velocity
p =    Coef_p * vT ;
u =    - (1.0 / bedParams.Kappa) * (Coef_p * vTz) ;  // pressure gradient determines velocity
u = Coef_u * vT ;  // makes my plotting easy  (this is just curve fitting)

// Nodal Values of Gas phase compostion, adsorbed/equil Concentrations in solid
// and adsoption rate into solid

for n in 1:Nc-1 loop              
      y[:, n] =   Coef_y[:, n] * vT    ;
end for ;

for i in 1:N loop
      // y[i, Nc] = 1 .- min(1, sum ( y[i, n] for n in 1 : Nc-1 )) ; 
      //y[i, Nc] = 1 .-  sum ( y[i, n] for n in 1 : Nc-1 ) ; 
       sum ( y[i, n] for n in 1 : Nc ) = 1 ; 
end for ; 

for n in 1:Nc loop              
      Q[:, n] =  Coef_Q[:, n] * vT   ;
    for m in 1:N loop
           Qeq[m, n] = max(0, 
                           bedParams.Qs[n] * ( bedParams.B[n] *  p[m] * y[m, n] ) 
                           / ( bedParams.Tb + p[m] * sum ( bedParams.B[j] *  y[m, j] for j in 1 : Nc) ));
          S[m, n]    =   bedParams.Km[n] * (Qeq[m, n] - max(0, Q[m, n])) ;  // Rate of adsorption
    end for ;
end for ;


// Component Balance Eqns.
 for  n in  1:Nc-1 loop       // nth Component Balance
     for m in 1:N-2 loop      // interior Collocation
         p[m] * ( der (Coef_y[:, n]) * vT[:, m]                          // p * (dy_i/dt
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
       der (Coef_p) * vT[:, m]                  // dp/dt  
            + u[m] * Coef_p * vTz[:, m]        // + u * dp/dz
            + p[m] * ( - (1.0 / bedParams.Kappa) * Coef_p * vTzz[:, m] )        // + p * du/dz
               + bedParams.Epsilon                                      // (1-e)/e
                    *  bedParams.Tb * sum ( S[m, j] for j in 1 : Nc )   //    T     * sum (dQ_i/dt)
               = 0 ; 
     end for;


// Adsorption rate into the solid
 for  n in  1:Nc loop       // nth Component 
   for m in 1:N loop      // interior Collocation
        der (Coef_Q[:, n]) * vT[:, m]                        // dQ_i/dt
                         - S[m, n]                           // - Km_i * (Qi_e - Qi)
          = 0 ;   
     end for;
 end for; // end of all component balances

// Boundary Conditions

 for  n in  1:Nc-1 loop       // nth Component 

    Coef_y[:, n] * vTz[:, N-1] =  bedParams.Pe * max(u[N-1], 0) * (y[N-1, n] - yin_in[n]) ; // bed inlet
    Coef_y[:, n] * vTz[:, N]   =  bedParams.Pe * min(u[N], 0) * (y[N, n] - yin_out[n]) ; // bed outlet 

 end for;

// Pressure Boundary Conditions
//             Coef_p[:] * vT[:, N-1] = p_in   ;    // Inlet  pressure 
//             Coef_p[:] * vT[:, N]   = p_out    ;    // outlet  pressure 

         0 = if (inlet.m_flow == 0.0) then Coef_p[:] * vTz[:, N-1]
             else Coef_p[:] * vT[:, N-1] - p_in ;    // Inlet  pressure 
          
          0 = if (outlet.m_flow == 0.0) then Coef_p[:] * vTz[:, N]
              else Coef_p[:] * vT[:, N] - p_out ;    // outlet  pressure 

//        Coef_p[:] * vTz[:, N]   = 0;    //  (zero gradient = no flow)

end partAdsorber;



/*
    Coef_y[:, n] * vTz[:, N-1] =  if(u[N-1] > 0) then // (p_in - p[1]), 
                                      bedParams.Pe * u[N-1] * (y[N-1, n] - yin_in[n])
                                  else  0 ; // bed inlet

    Coef_y[:, n] * vTz[:, N]   =  if(u[N] > 0) then   // (p[N] - p_out), 
                                    0 
                                  else bedParams.Pe * u[N] * (y[N, n] - yin_out[n]) ; // bed outlet 

    Coef_y[:, n] * vTz[:, N-1] =  regStep(u[N-1]-uTol, // (p_in - p[1]), 
                                      bedParams.Pe * u[N-1] * (y[N-1, n] - yin_in[n]),  0, uTol) ; // bed inlet

    Coef_y[:, n] * vTz[:, N]   =  regStep(u[N]+uTol,    // (p[N] - p_out), 
                                    0,  bedParams.Pe * u[N] * (y[N, n] - yin_out[n]), uTol) ; // bed outlet 

    Coef_y[:, n] * vTz[:, N-1] =  ThermoS.Math.regStep(u[N-1]-uTol, // (p_in - p[1]), 
                                      bedParams.Pe * u[N-1] * (y[N-1, n] - yin_in[n]),  0, uTol) ; // bed inlet

    Coef_y[:, n] * vTz[:, N]   =  ThermoS.Math.regStep(u[N]+uTol,    // (p[N] - p_out), 
                                    0,  bedParams.Pe * u[N] * (y[N, n] - yin_out[n]), uTol) ; // bed outlet 
             Coef_p[:] * vT[:, N-1] = p_in   ;    // Inlet  pressure 
             Coef_p[:] * vT[:, N]   = p_out    ;    // outlet  pressure 

    Coef_y[:, n] * vTz[:, N-1] =  if(u[N-1] > 0) then 
                                      bedParams.Pe * u[N-1] * (y[N-1, n] - yin_in[n])
                                  else  0 ; // bed inlet

    Coef_y[:, n] * vTz[:, N]   =  if(u[N] >= 0) then   
                                    0 
                                  else bedParams.Pe * u[N] * (y[N, n] - yin_out[n]) ; // bed outlet 

// Smoothedout B.C changes during flow reversal

    Coef_y[:, n] * vTz[:, N-1] =  ThermoS.Math.regStep(p[N-1] - p[1], 
                                          bedParams.Pe * u[N-1] * (y[N-1, n] - yin_in[n]), 0.0, pTol) ; // bed inlet

    Coef_y[:, n] * vTz[:, N]   =  ThermoS.Math.regStep(p[N] - p[N-2],
                                          bedParams.Pe * u[N] * (y[N, n] - yin_out[n]), 0.0, pTol) ; // bed outlet 
*/
