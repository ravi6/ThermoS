within Rig;
model Tee
/*
*  A Three Port Mixer  with finite volume
*/

    replaceable package Medium = PartialMixtureMedium ;

//  Parameters
  parameter    Volume   	   vol   	= 1e-4 ;   // Tank Volume (m3)
  parameter    Boolean             Adiabatic  	= true ;   // default is adiabatic 
  parameter    Medium.Temperature  Tset 	= 300   ;   //  at 300K
  parameter    Integer		   N 		= 3     ;   // Number of Ports

  FluidPort port[N] (redeclare each package Medium = Medium) ;

    Medium.Temperature		T      ;
    Medium.AbsolutePressure	p      ;
    HeatFlowRate		Q_in   ;
    Medium.BaseProperties       medium ;

    Mass			M (start=1e-4)	;	// Mass of Gas in the vessal 
    InternalEnergy              U               ;       // Enthalpy of contents

  equation

    if (Adiabatic) then
      Q_in = 0 ;
    else
      medium.T = Tset ;
    end if;


     // Mass and Component Balance
     M = vol * medium.d ;
     der(M) = sum (port[:].m_flow) ;  // Mass Balance

     // Enthalpy Balance
     U = M * medium.u ;
     der (U) =  sum (port[i].m_flow * actualStream(port[i].h_outflow)
                                  for i in 1:N)
    	                 + Q_in; 

     // Tee Junction well mixed
      for i in 1:N loop
        port[i].Xi_outflow = medium.Xi ; 
	port[i].h_outflow  = medium.h  ;
        port[i].p = medium.p ;
      end for;

      T = medium.T ;
      p = medium.p ;  
      medium.Xi = Medium.reference_X [1:Medium.nXi] ;

end Tee ;
