model plant
/*
  Author: Ravi Saripalli
  	 12th May 2014 
*/
  import ThermoS.Uops.*;
  import ThermoS.Types.*;
  import ThermoS.Media.MyGas;

  constant Real AirComp[2] = {0.767,0.233}; // reduced Xi
  parameter Integer Nseg =20  ;
  parameter Boolean CoCurrent = false ;

  Reservoir     res1	(redeclare package Medium = MyGas, p = 5e5, T = 300, Xi = AirComp); // Reservoir 1
  Reservoir 	res2	(redeclare package Medium = MyGas, p = 5e5, T = 600, Xi = AirComp); // Reservoir 2
  Reservoir 	res3	(redeclare package Medium = MyGas, p = 1e5, T = 300, Xi = AirComp); // Reservoir 2
  Reservoir 	res4	(redeclare package Medium = MyGas, p = 1e5, T = 300, Xi = AirComp); // Reservoir 2
  HxMultiSeg 	hx	(redeclare package Medium_c = MyGas, redeclare package Medium_h = MyGas,
                                           nseg = Nseg, Awf_h = 10/Nseg, Awf_c = 10/Nseg,
                                                     holdup_h = 50.0/Nseg, holdup_c = 50.0/Nseg,
                                                     w_m = 1.0/Nseg, 
                                                     cf_h = 1e-3*sqrt(Nseg), cf_c = 1e-3*sqrt(Nseg)); 

equation

     connect (res1.port, hx.portA_c);
     connect (res3.port, hx.portB_c);
     if CoCurrent then
       connect (res2.port, hx.portA_h);
       connect (res4.port, hx.portB_h);
     else
       connect (res4.port, hx.portA_h);
       connect (res2.port, hx.portB_h);
     end if;        
    


initial equation
   for i in 1:Nseg loop
     hx.seg_h[i].Tf = 300 ;
     hx.seg_h[i].Tf = hx.seg_c[i].Tf;
     hx.Tw[i] = 300 ;
   end for;

end plant;
