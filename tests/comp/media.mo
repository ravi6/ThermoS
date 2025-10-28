  model  plant
  import ThermoS.Types.* ;
  import ThermoS.Media.MyGas ;
  replaceable package  medium = MyGas ;

    // Output variables for medium properties
    Real h,s,T, x(start=0, fixed=true), mflow ;
    medium.BaseProperties  bp, bp2 ;

  equation
    bp.p = 1e5 ;
    bp.T = 300 ;
    bp.Xi = {0.7, 0.2} ;
    h = medium.specificEnthalpy (bp.state);
    s = medium.specificEntropy (bp.state) ;

//    bp2.p = 4e5 ;
    mflow = 0.001 ;
    bp2.Xi = bp.Xi ;
    bp2.h = medium.specificEnthalpy (medium.setState_psX(bp2.p, s, bp2.X)) ;
    bp2.T = T ;
    x = - mflow * (h - bp2.h) ;
    der(x) = 30 ;
     assert(x>0, "T=" + String(T), level = AssertionLevel.warning);
    
  end plant ;
