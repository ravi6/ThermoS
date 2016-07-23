within ThermoS.Uops.Valves;

model Valve                        // A control Valve
    extends ThermoS.Uops.Valves.partialValve  ;

parameter Vchar vchar = Vchar.Linear             ;  // Valve Charachteristics (Linear by default)

Percent      po    (start=50)            ;               // Valve % Open  
Fraction     charF (start=1.0)           ;               // Characteristic Multiplier


equation

  if (vchar == Vchar.Linear) then
     charF = po / 100 ;
  elseif (vchar == Vchar.FastActing) then 
     charF = (po/100)^0.5 ; 
  elseif (vchar == Vchar.EquiPercent) then
     charF =  35^(po/100-1) ; // base could be varied from 20 to 50
  end if;

    //inlet.m_flow =  cv * charF  * sqrt(max(0,med.d)) * noEvent(regRoot(inlet.p - outlet.p, dpTol )) ;
    inlet.m_flow =  cv * charF  *  sqrt(med.d) *  regRoot(inlet.p - outlet.p, dpTol ) ;

end Valve;









  /* Note the regRoot function regularizes sign(x) * sqrt(abs(x)) function near at x=0
       using default delta 0.01. ie. when pressure drop is less than 1% of inlet pressure
       the deviation of mass flow rate  from square root law is about 0.25%
     Homotopy is used to assit convergence for very small pressure drop .. ie. linear model
     to non-linear model morphing.
  */


//    assert(med.d>0,  "Screwed up Compo" 
//                                     +  ThermoS.Util.strVec(med.state.X) 
//                                     , AssertionLevel.warning);

 //   inlet.m_flow = cv * charF  * sqrt( max(0,Medium.density(state)) )
 //                               *   regRoot(inlet.p - outlet.p, dpTol ) ;

   // bypassing Medium.density algorithmeic call that is giving some angst to the nonlinear solver 
//    inlet.m_flow = cv * charF  * sqrt( max(0, (state.p/((Medium.data.R * state.X) * state.T))) )
//                                *   regRoot(inlet.p - outlet.p, dpTol ) ;
/*
                      * homotopy( (1 - outlet.p / inlet.p) , 
                                  regRoot(1 - outlet.p / inlet.p)
                                );  // dumping due index reduction issues with adsorber
*/

