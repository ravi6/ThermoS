// Version 1.0
// Date:   9th December 2013
within ThermoS.Uops.Tanks;
model portMixer
/*
*  A Three Port Mixer  with finite volume
*/

    replaceable package Medium = PartialMixtureMedium ;

//  Parameters
  parameter    Volume   	   vol   	= 1e-6 ;   // Tank Volume (m3)
  parameter    Boolean             Adiabatic  	= false ;   // default is Isothermal
  parameter    Medium.Temperature  Tset 	= 300   ;   //  at 300K
  parameter    Integer		   N 		= 2     ;   // Number of Ports

  FluidPort port[N] (redeclare each package Medium = Medium) ;

    Medium.Temperature		T      ;
    Medium.AbsolutePressure	p      ;
    HeatFlowRate		Q_in   ;
    Medium.BaseProperties       medium ;

    MassFraction                Xi[Medium.nXi]  ;
    Mass			M		;	// Mass of Gas in the vessal 
    Enthalpy                    H               ;       // Enthalpy of contents

  equation

    if (Adiabatic) then
      Q_in = 0 ;
    else
      medium.T = Tset ;
    end if;


     // Mass and Component Balance
     M = vol * medium.d ;
     der(M) = sum (port[:].m_flow) ;  // Mass Balance

     // Component Balance
     for j in 1:Medium.nXi loop
         der (M * Xi[j]) = 
              sum (actualStream(port[i].Xi_outflow[j]) * port[i].m_flow
                                   for i in 1:N) ;
     end for;

     // Enthalpy Balance
     H = vol * medium.d * medium.h ;
     der (H) =  sum (port[i].m_flow * actualStream(port[i].h_outflow)
                                  for i in 1:N)
    	               + vol * der(medium.p)  + Q_in; 

     // Assume gas in tank is well mixed 
      for i in 1:N loop
        port[i].Xi_outflow = medium.Xi ; 
	port[i].h_outflow  = medium.h  ;
        port[i].p = medium.p ;
      end for;

      T = medium.T ;
      p = medium.p ;  
      Xi = medium.Xi ;

end portMixer;
