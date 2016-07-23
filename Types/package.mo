//Author: Ravi Saripalli
within ThermoS;
  package Types	"Contains all ThermoS Typed Variables"

/* we can dispense with these now that we have SIunits 
 import Modelica.SIunits.*;
type Pressure    = Real(unit="Pa")    ;
type Temperature = Real(unit="C")     ;
type MassFlow    = Real(unit="kg/s")  ;
type EnthFlow    = Real(unit="J/s")   ;
type Work        = Real(unit="W")     ;
type Length      = Real(unit="m")     ;
type Density     = Real(unit="kg/m3") ;
type Area        = Real(unit="m2")    ;
type MassFlow = MassFlowRate;    // This is how I can make my own types
                                //  from SIunits Namespace
				//  Just an illustraton
*/

type Percent = Real(unit="%", min=0 , max=100) ;
type Fraction = Real(min=0, max=1.0) ;

end Types;
