//Author: Ravi Saripalli
within ThermoS.Uops.PSA;

function moleToMass "Mole fractions to Mass fractions"

    input MoleFraction      moleFrac[:]   "Mole fractions of mixture";
    input MolarMass         Mw[:]         "Molar masses of components (all components)";
    output MassFraction     massFrac[size(moleFrac,1)]   "Mass fractions";


  protected
     Integer Nx = size(moleFrac,1)  ;
     Integer Nc = size(Mw,1);
     MolarMass    Mmix  ;

  algorithm

    Mmix := 0 ;
    for i in 1:Nx loop
       Mmix :=  Mmix + moleFrac[i] * Mw[i]  ;
    end for;

    if ( Nx < Nc ) then 
       Mmix := Mmix + (1 - sum(moleFrac)) * Mw[Nc] ;
    end if;

    for i in 1:Nx loop
      massFrac[i] := moleFrac[i] * Mw[i] / Mmix;
    end for;

end moleToMass;
