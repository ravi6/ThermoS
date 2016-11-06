within ThermoS.Media;

package MyAir   "Specific ideal Gas Mixture of N2 and O2"

import Modelica.Media.IdealGases.Common.MixtureGasNasa;
import Modelica.Media.IdealGases.Common.SingleGasesData ;
import Modelica.Media.IdealGases.Common.FluidData;

extends MixtureGasNasa (
	 data = {SingleGasesData.N2, // note data is of type  DataRecord[:]
		 SingleGasesData.O2},
	 fluidConstants = {FluidData.N2,
			   FluidData.O2}, 
	 substanceNames = {"Nitrogen", "Oxygen"},
         reducedX = true, 
         reference_X = {0.768, 0.232},

         //preferredMediumStates = true, // p,T, Xi are the preferred states I think
         // The above line gives trouble with CODE GEN why ... especially with my singleBedModel "Case3"
         // Prepare a bug report ...
// Modifying start values and Ranges 
         Density(start=1, nominal=1, min=1e-3, max=10),
         AbsolutePressure(start=1e5, min=1e3, max=50e5, nominal=1e5),
         Temperature(start=300, min=200, max=2000, nominal=300),
         ThermodynamicState(p(start=1e5), T(start=300), X(start=reference_X)),
         MassFraction(each start=0.5),
         SpecificEnthalpy(start=300e3, nominal=300e3) 
); //end of Extending MixtureGasNasa
         constant MoleFraction reference_Xm[2] = {0.79, 0.21} ;  // Molar fractions
end MyAir;
