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
         reference_X = {0.7, 0.2, 0.1}
);
end MyGas; 
