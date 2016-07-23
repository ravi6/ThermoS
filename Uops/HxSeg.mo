within ThermoS.Uops;

/*  A generic heat exchanger Segment 
      Exchanges heat between two fluids 
      separated by a wall. Heat loss to other surrounds
      from each fluid can be specified.
*/

model HxSeg
  replaceable package Medium_h  =  PartialMixtureMedium ;
  replaceable package Medium_c  =  PartialMixtureMedium ;

  FluidPort 	portA_h (redeclare package Medium = Medium_h)  ; 
  FluidPort 	portB_h (redeclare package Medium = Medium_h)  ; 
  FluidPort 	portA_c (redeclare package Medium = Medium_c)  ; 
  FluidPort 	portB_c (redeclare package Medium = Medium_c)  ; 

  parameter Real	cf_h	= 0.001     ; // Pressure Loss Coeff. (m^3/Pa^0.5))     
  parameter Real	cf_c	= 0.001     ; // Pressure Loss Coeff. (m^3/Pa^0.5))     
  parameter Area	Awf_h	= 1.0	  ; // Wall to fluid heat transfer area
  parameter Area	Awf_c	= 1.0	  ; // Wall to fluid heat transfer area
  parameter CoefficientOfHeatTransfer	hwf_c	= 150	; // Wall to fluid heat transfer coeff. 
  parameter CoefficientOfHeatTransfer	hwf_h	= 150	; // Wall to fluid heat transfer coeff. 
  parameter Mass		        w_m	= 1.0	; // Mass of heat transfer walls
  parameter SpecificHeatCapacity        w_cp    = 420	; // Specific heat of wall material
  parameter Volume			holdup_h  = 50  ; // Heater fluid holdup 
  parameter Volume			holdup_c  = 50  ; // Heater fluid holdup 

  Temperature			Tw		    ; // Wall temperature (K)
  Heat				Qwf_h		    ; // Heat tranfer from wall to fluid
  Heat				Qwf_c		    ; // Heat tranfer from wall to fluid
  Heat				Qloss_c = 0	    ; // Heat loss to external environment from cold side 
  Heat				Qloss_h = 0	    ; // Heat loss to external environment from hot side 
  HxFluidSeg   seg_h(redeclare package Medium = Medium_h, cf = cf_h, holdup = holdup_h);
  HxFluidSeg   seg_c(redeclare package Medium = Medium_c, cf = cf_c, holdup = holdup_c);

  equation
      connect   (seg_h.portA,  portA_h)  ;
      connect   (seg_h.portB,  portB_h)  ;
      connect   (seg_c.portA,  portA_c)  ;
      connect   (seg_c.portB,  portB_c)  ;
	Qwf_h = hwf_h * Awf_h * (Tw - seg_h.Tf)  ;  // wall to hot fluid
	Qwf_c = hwf_c * Awf_c * (Tw - seg_c.Tf)  ;  // wall to cold fluid
        w_m * w_cp * der(Tw) =  - (Qwf_h + Qwf_c);  //  Wall thermal inertia
                                     /* ignoring wall conduction resistance */
       /* Heat loss to the external world */
	seg_c.Qin = Qwf_c - Qloss_c ;
        seg_h.Qin = Qwf_h - Qloss_h ;
end HxSeg;
