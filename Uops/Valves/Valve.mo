within ThermoS.Uops.Valves;

model Valve                        // A control Valve
    extends ThermoS.Uops.Valves.partialValve  ;

// import ThermoS.Math.sMax0 ;

parameter Vchar vchar = Vchar.Linear             ;  // Valve Charachteristics (Linear by default)
parameter Fraction pratChoke  = 0.5              ;  // Maximum pressure ratio at choking point
parameter Boolean  Compressible = true           ;  // Default to compressible flow

Percent      po    (start=0.0)           ;               // Valve % Open  
Fraction     charF (start=1.0)           ;               // Characteristic Multiplier
Fraction     prat  (start=1.0)           ;

equation

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
          inlet.m_flow = noEvent (
                           if ( prat > 1  or prat < 1) then
                             cv * max(0,charF) *  sqrt(max(0,med.d)) 
                                * sqrt(max(inlet.p, outlet.p))
                                * sign(inlet.p - outlet.p)  
                                * regRoot(1 - max(prat, 0.5), dpTol )
                            else 0 ) ; 
    else
       inlet.m_flow = cv * charF  *  sqrt(med.d * inlet.p)
                        * regRoot2(1 - outlet.p/inlet.p, dpTol) ; 
    end if; 

    
end Valve;
