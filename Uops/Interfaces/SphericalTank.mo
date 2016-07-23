within ThermoS.Uops.Interfaces;

model SphericalTank

   import    Modelica.Constants.* ;

   parameter Length  d   = 1.0 ;  // Tank  Diameter (m)
       
   Volume  fvol   ;  // Liquid volume in the tank (m3)
   Volume  tvol   ;  // Total tank volume (m3)
   Length  level  ;  // Tank Level (m)
   Percent pFull  ; // %Tank full
   Area    wettedArea ;  // Area wetted by the fluid
   Area    interfaceArea ;  // gas-liquid interfacial area

    equation
             fvol =  pi * level * level * (1.5 * d  - level) / 3 ;
             tvol = pi * d * d * d / 6     ;
             pFull =  100 * fvol / tvol ;
             wettedArea =    pi * d * level   ;
             interfaceArea = pi * ( d * d / 4 
                                      + (level - d/2)^2 ) ;
end SphericalTank ;
