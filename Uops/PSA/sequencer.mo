within ThermoS.Uops.PSA;
model sequencer
/* This implements discrete events
      to simulate Adsorber Valving Sequence
*/ 
import  Modelica.SIunits.Time ;

parameter Time pTime = 5 ;    // pressurization Time (sec)
parameter Time fTime = 10 ;    //  feeding  Time (sec)
parameter Time puTime = 5 ;    // purging Time (sec)
parameter Time sTime = 1  ;    // Descrete Clock Sample Time (sec)


Boolean pressurizing, feeding, purging ;
Boolean pressurize, feed,  purge;
Integer pcount, fcount, pucount ;

Integer tics  ;
Boolean clocktic  ;

 equation

// Digital Clock
    clocktic = sample(0,sTime) ;  
    when clocktic then 
       tics = pre(tics) + 1 ;
    end when;
  
 // State Transition Logic     
  when pressurizing then
    pcount = tics ;     // reset pressureization counter
  end when;
  
  when feeding then
    fcount = tics ;
  end when;
  
  when purging then 
    pucount = tics ;
  end when ;
  
 // Time based transitions
    
when clocktic  then
   feed = (pre(pressurizing) and ((tics - pre(pcount)) * sTime >= pTime)) ;
   purge = (pre(feeding) and ((tics - pre(fcount)) * sTime>= fTime)) ;
   pressurize = (pre(purging) and ((tics - pre(pucount)) * sTime >= puTime)) ;
end when ;

when {feed, pressurize, purge} then
   if (pressurize and pre(purging) ) then
     pressurizing = true ; feeding = false ; purging = false ; 
   elseif (feed and pre(pressurizing)) then
     pressurizing = false ; feeding = true ; purging = false ; 
   elseif (purge and pre(feeding) ) then
     pressurizing = false ; feeding = false ; purging = true ; 
   else
     pressurizing = pre(pressurizing) ; 
     feeding = pre(feeding) ; purging = pre(purging) ; 
   end if;
end when;
 
  initial equation
  pcount = 0 ; pucount = 0 ; fcount = 0 ; tics = 0 ;
//  pressurizing = false ; feeding = false ; purging = true;
  pressurize = false ; feed = false ; purge = false ;
  
end sequencer;
