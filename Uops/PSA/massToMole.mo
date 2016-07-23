within ThermoS.Uops.PSA;
//Author: Ravi Saripalli

// Some interface functions I need in this package (hopefully not may)
//   These functions really should be unnecssary if only we have a more generic
// functions moleToMassFractions and massToMoleFractions functions that can
//   handle reduced index composition formulations (with NS-1 compositions specified)

function massToMole "Mass fractions to Mole fractions"
// This function is a bit better than that is provided in the
//    partialMixture package. We can use it does not rely on complete
//    composition array.   Will handle structurally independent composition.
//    But of course you still need to   give full Molec.Wt. array 

   input MassFraction     massFrac[:] "Mole fractions of mixture"; 
   input MolarMass        Mw[:]       "Molar masses of components (all components)";
   output MoleFraction    moleFrac[size(massFrac,1)] "Mole fractions"; 

   protected
     Integer Nx = size(massFrac, 1) ; 
     Integer Nc = size(Mw, 1)       ;
     MolarMass    Moles[size(Mw,1)]     ;  

  algorithm

        for i in 1 : Nx loop
           Moles[i] := massFrac[i] / Mw[i] ; 
        end for;

       if ( Nx < Nc ) then
         Moles[Nc] := (1 - sum(massFrac)) / Mw[Nc] ;
       end if;

      for i in 1:Nx loop
         moleFrac[i] := Moles[i]/sum(Moles);
      end for;

end massToMole;
