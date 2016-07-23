within ThermoS.Uops.PSA;
//Author: Ravi Saripalli

record bedParamsRec

  constant Integer  Nc = 2      ; // Number species in Gas media 
  constant Real R = 8314.5      ; // Universal Gas Constant Pa. m3/(C kgmol)

  Real Pref = 1e5               ; // Reference pressure (Pa)
  Real Tref = 300               ; // Reference Temperature (K) 
  Real Uref = 1.0               ; // Refernece Gas Velocity (m/s)
  
  Real Diff = 1e-3              ; // Axial dispersion Coeff  (m2/sec) 
  Real mu   = 2.0e-5            ; // Air viscosity (kg/m.s)
  Real dp   = 1.0e-3            ; // Adsorbent particle size (m)
  Real Tbed = 300               ; // Adsorber temperature (K)

  Real voidage = 0.5            ; // Bed voidage
  Real L    = 0.25              ; // Adsorber Bed length (m)
  Real dia  = 0.2               ; // Adsorber Bed diameter(m)
  
  Real qs[Nc]  = {5.26, 5.26}        ; // Langmuir Saturation Concentration (kgmol/m3)
  Real hk[Nc]  = {14.8, 4.7}         ; // Henry's Constants (dim. less) (= qs * b)
  Real k_m[Nc] = {19.7, 62.0}        ; // External Mass transfer rate s (1/s)
  
  // Derived quantities
  Real csArea = 3.147*dia*dia/4 ; // Bed CrossSectional Area (m2)   

  Real Cref    = Pref / (R * Tref)             ; // Reference Molar Concentration (kgmol/m3)
  Real b[Nc]   = hk ./ qs                      ; // Units m3/kgmol
  Real Kozney  = 180 * mu * (1 - voidage)^2    
                          / ( dp * dp * voidage^3 )    ; // KozneyCarman Eqn. Coeff. (pressure grad/vel)
  
  // These are the basic params the core model needs
  Real Epsilon = (1 - voidage) / voidage     ;
  Real Km[Nc]  = k_m * (L / Uref)            ;// Dim.less Mass trasnfer rate 
  Real Kappa   = Kozney * L * Uref / Pref    ;// dp*/dz* = - Kappa* v*   (dimensionless)
  Real Pe      = L * Uref / Diff             ;// Peclet number
  Real B[Nc]   = b * Cref                    ;// Normalized b
  Real Qs[Nc]  = qs / Cref                   ;// Normalized Saturation Concentration 
  Real Tb      = Tbed / Tref                 ;// Adsorber temperature (isothermal)

end bedParamsRec;
