within ThermoS.Media;

package MyGas   "Specific ideal Gas Mixture"

import Modelica.Media.IdealGases.Common.MixtureGasNasa;
import Modelica.Media.IdealGases.Common.SingleGasesData ;
import Modelica.Media.IdealGases.Common.FluidData;

extends MixtureGasNasa (
	 data = {SingleGasesData.N2, // note data is of type  DataRecord[:]
		 SingleGasesData.O2,
		 SingleGasesData.CO2},
	 fluidConstants = {FluidData.N2,
			   FluidData.O2, 
			   FluidData.CO2}, 
	 substanceNames = {"Nitrogen", "Oxygen", "CarbonDioxide"},
         reducedX = true,
         reference_X = {0.7, 0.2, 0.1},

// Modifying start values and Ranges 
         Density(start=1, nominal=1),
         AbsolutePressure(start=1e5, min=1e3, max=50e5, nominal=1e5),
         Temperature(start=300, min=200, max=2000, nominal=300) //,
         //MassFraction(start=0.333333), MoleFraction(start=0.333333)
);
end MyGas; 
