//Author: Ravi Saripalli
within ThermoS;
  package Types	"Contains all ThermoS Typed Variables"

type Percent = Real(unit="%", min=0 , max=100) ;
type Fraction = Real(min=0, max=1.0) ;

end Types;
