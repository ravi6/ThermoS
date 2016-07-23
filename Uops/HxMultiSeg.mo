within ThermoS.Uops;
/*  A generic heat exchanger 
      Exchanges heat between two fluids 
     with interconnected HxSegments
*/

model HxMultiSeg
  replaceable package Medium_h  =  PartialMixtureMedium ;
  replaceable package Medium_c  =  PartialMixtureMedium ;

  FluidPort 	portA_h (redeclare package Medium = Medium_h)  ; 
  FluidPort 	portB_h (redeclare package Medium = Medium_h)  ; 
  FluidPort 	portA_c (redeclare package Medium = Medium_c)  ; 
  FluidPort 	portB_c (redeclare package Medium = Medium_c)  ; 

  parameter Integer     nseg=1	;  // Number of segments
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
  parameter Heat 			Qloss_c = 0	    ; // Heat loss to external environment from cold side 
  parameter Heat 			Qloss_h = 0	    ; // Heat loss to external environment from hot side 

// Things that vary along the segment
  Temperature[nseg]		Tw		    ; // Wall temperature (K)
  Heat[nseg]			Qwf_h		    ; // Heat tranfer from wall to fluid
  Heat[nseg]		        Qwf_c		    ; // Heat tranfer from wall to fluid
  HxFluidSeg[nseg]    seg_h( redeclare each package Medium = Medium_h, 
                           each cf = cf_h,  each holdup = holdup_h);
  HxFluidSeg[nseg]  seg_c( redeclare each package  Medium = Medium_c, 
                          each cf = cf_c, each holdup = holdup_c);
  equation

// Inter connect segments in line (a--b a--b a--b etc.)
     for i in  1:nseg-1 loop
      connect   (seg_h[i].portB,  seg_h[i+1].portA)  ;
      connect   (seg_c[i].portB,  seg_c[i+1].portA)  ;
     end for;

// Hook up end ports of extreme segments
      connect   (seg_h[1].portA,     portA_h)  ;
      connect   (seg_h[nseg].portB,  portB_h)  ;
      connect   (seg_c[1].portA,     portA_c)  ;
      connect   (seg_c[nseg].portB,  portB_c)  ;

// Segment-wise equations
     for i in  1:nseg loop
	Qwf_h[i] = hwf_h * Awf_h * (Tw[i] - seg_h[i].Tf)  ;  // wall to hot fluid
	Qwf_c[i] = hwf_c * Awf_c * (Tw[i] - seg_c[i].Tf)  ;  // wall to cold fluid
        w_m * w_cp * der(Tw[i]) =  - (Qwf_h[i] + Qwf_c[i]);  //  Wall thermal inertia
                                     /* ignoring wall conduction resistance */
       /* Heat loss to the external world */
	seg_c[i].Qin = Qwf_c[i] - Qloss_c ;
        seg_h[i].Qin = Qwf_h[i] - Qloss_h ;
     end for;
end HxMultiSeg;
