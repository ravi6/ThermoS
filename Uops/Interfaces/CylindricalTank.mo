within ThermoS.Uops.Interfaces;

model CylindricalTank

   import    Modelica.Constants.* ;

   parameter Length  h   = 1.0 ;  // Tank Height (m)
   parameter Length  d   = 1.0 ;  // Tank Diameter (m)
       
   Volume  fvol   ;  // Liquid volume in the tank (m3)
   Volume  tvol   ;  // Total tank volume (m3)
   Length  level  ;  // Tank Level (m)
   Percent pFull  ; // %Tank full
   Area    wettedArea ;  // Area wetted by the fluid
   Area    interfaceArea ;  // gas-liquid interfacial area

    equation
             fvol = pi * d * d  * level / 4 ;
             tvol = pi * d * d  * h / 4     ;
             pFull =  100 * fvol / tvol ;

             wettedArea =  ( pi * d * d / 4 ) + ( pi * d * level ) ;
             interfaceArea = pi * d * d / 4 ;
        
end CylindricalTank ;
