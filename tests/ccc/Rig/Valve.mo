within Rig;

model Valve                        // A control Valve

  replaceable package Medium = PartialMixtureMedium;
  FluidPort inlet(redeclare package Medium = Medium);
  FluidPort outlet(redeclare package Medium = Medium);

  type Vchar = enumeration ( Linear         "Linear Valve",
                             FastActing     "Fast Acting Valve",
                             EquiPercent    "Equi-Percent Valve"
                           ) "Enumeration Defining Valve Behaviour" ; 

  parameter Vchar vchar = Vchar.Linear  ;  // Valve Charachteristics (Linear by default)
  parameter Fraction pratChoke  = 0.5   ;  // Maximum pressure ratio at choking point
  parameter Boolean  Compressible = true  ;  // Default to compressible flow

  // Valve Coefficient   (1 m3/s @ one bar differential with 1kg/m3 density)
  parameter Real cv=1.0/sqrt(1e5);
  parameter Real dpTol = 100 ;  // Pressure Drop Tolerance used for state and flow transition

  Percent      po    (start=0.0)  ; // Valve % Open  
  Fraction     charF (start=1.0)  ; // Characteristic Multiplier
  Fraction     prat  (start=1.0)  ; // Pressure Ratio
  Fraction     Y     (start=1.0)  ; // Compressibility Factor

  Medium.BaseProperties medium  ;  
 
  equation

  // Mass balance
  inlet.m_flow + outlet.m_flow = 0;

  outlet.h_outflow  = inStream(inlet.h_outflow); // flow is  inlet to outlet
  inlet.h_outflow   = inStream(outlet.h_outflow); // reverse flow case

  medium.p  = max(inlet.p, outlet.p);
  medium.h  = inlet.h_outflow  ;  
  medium.Xi = Medium.reference_X [1:Medium.nXi] ;


  if (vchar == Vchar.Linear) then
     charF =  po / 100 ;
  elseif (vchar == Vchar.FastActing) then 
     charF = (po/100)^0.5 ; 
  elseif (vchar == Vchar.EquiPercent) then
     charF =  35^(po/100-1) ; // base could be varied from 20 to 50
  end if;

    prat =  min(inlet.p, outlet.p) / max(inlet.p, outlet.p) ;

    // Make Valve equaiton  continuous and differentiable both near zero flows
    if (Compressible) then
          Y = 1 - (1 - prat) / ( 3 * (1 - pratChoke) ) ;
          inlet.m_flow = noEvent (
                           if ( prat > 1  or prat < 1) then
                             cv * Y * max(0,charF) *  sqrt(max(0,medium.d)) 
                                * sqrt(max(inlet.p, outlet.p))
                                * sign(inlet.p - outlet.p)  
                                * regRoot(1 - max(prat, 0.5), dpTol )
                            else 0 ) ; 
    else
       inlet.m_flow = cv * charF  *  sqrt(medium.d * inlet.p)
                        * regRoot2(1 - outlet.p/inlet.p, dpTol) ; 
       Y = 1  ; // redundant but needed in new version
    end if; 
    
    inlet.Xi_outflow = medium.Xi ;
    outlet.Xi_outflow = medium.Xi ;
end Valve;
