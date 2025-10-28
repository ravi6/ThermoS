  model  plant
  import ThermoS.Types.* ;
  import ThermoS.Media.MyGas ;
  replaceable package  medium = MyGas ;

    // Output variables for medium properties
    Real h ,T;
    medium.BaseProperties  bp ;

  equation
    bp.p = 1e5 ;
    bp.T = 300 ;
    bp.X = {0.7, 0.2,0.1} ;
    h = medium.specificEnthalpy (bp.state) ;
  end plant ;
