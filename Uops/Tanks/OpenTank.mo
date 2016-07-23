within ThermoS.Uops.Tanks;

model OpenTank   " An Open Tank "
   extends ThermoS.Uops.Interfaces.CylindricalTank ; 
   extends SimpleTank ;

/* It has two connections
     inlet ... can be located at any level in the tank
     outlet ... fixed at bottom for now
     implements 
             Heat transfer into fluid through tank walls(wetted)
             Accounts for static head effects on inlet and outlet connectors
*/

//  Parameters
  parameter    Length	in_pos   = h   ; // at the top of the tank  
//  parameter    Length	out_pos  = 0   ;  // fixed for now  

  CoefficientOfHeatTransfer	hcoef   	    ; // Ambient to fluid heat transfer coeff. 
  Area	                        A_h	            ; //  heat transfer area
  Temperature			Ta		    ; // Ambient temperature (K)

  equation     

    outlet.p = Pa + rho * 9.8 * level  ; 
    assert(level <= h or inlet.m_flow > 0,  "Tank overflowing  at t="+String(time),
                      level = AssertionLevel.warning);
    assert(level > 0 ,  "Tank Empty  at t="+String(time),
                      level = AssertionLevel.error);

/* 
This sort of working but with 
    unrealistic wiggles near discontinuity 

   if (level > in_pos) then    // Inlet is below the tank level
   	inlet.p = Pa +  rho * 9.8 * (level - in_pos) ; 
        // asserts when the first argument is false
	assert(false, "Tank Level is above the Inlet at t="+String(time),
                      level = AssertionLevel.warning);
   else    // Inlet is below the tank level
   	inlet.p = Pa ;
	assert(false, "Tank Level is below the Inlet at t="+String(time),
                      level = AssertionLevel.warning);
   end if;

*/

// This is the correct method (first argument is the degree of smoothing)
     inlet.p = smooth(1,  if level > in_pos then 
   	                      Pa +  rho * 9.8 * (level - in_pos)  
                          else
                              Pa ) ;
  
	assert(level < in_pos,  "Tank Level is above the Inlet at t="+String(time),
                      level = AssertionLevel.warning);
	assert(level > in_pos, "Tank Level is below the Inlet at t="+String(time),
                      level = AssertionLevel.warning);
	assert((level<in_pos) or (level>in_pos), "Tank Level is at the Inlet at t="+String(time),
                      level = AssertionLevel.warning); // funny I did not hit this condition so far
       
        A_h = wettedArea ;
	Q_in = hcoef * A_h * (Ta - Tf)  ;

end OpenTank;
