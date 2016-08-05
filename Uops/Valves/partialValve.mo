within ThermoS.Uops.Valves;
partial model partialValve
import Modelica.Media.IdealGases.Common.MixtureGasNasa.h_TX;

  // A simple Valve
  replaceable package Medium = PartialMixtureMedium;

  // Specify that our Medium is used in/outlet
  FluidPort inlet(redeclare package Medium = Medium);
  FluidPort outlet(redeclare package Medium = Medium);

  // Valve Coefficient   (1 m3/s @ one bar differential with 1kg/m3 density)
  parameter Real cv=1.0/sqrt(1e5);
  parameter Real dpTol = 100 ;  // Pressure Drop Tolerance used for state and flow transition

 // Medium.ThermodynamicState state (T(start=300), p(start=1e5), X(start=Medium.reference_X));
  Medium.BaseProperties med; //(preferredMediumStates=true) ;  
 
equation

/* No change in Enthalpy in the Valve
    * inStream returns value when flow is into the device
    * Looks trivial for one to one connections. But is designed to handle
    * large number of connection branches without singularity
*/

  outlet.h_outflow  = inStream(inlet.h_outflow); // flow is  inlet to outlet
  inlet.h_outflow   = inStream(outlet.h_outflow); // reverse flow case

  inlet.Xi_outflow  = inStream(outlet.Xi_outflow); // pass composition un_altered
  outlet.Xi_outflow = inStream(inlet.Xi_outflow); // pass composition un_altered

  // Mass balance
  inlet.m_flow + outlet.m_flow = 0;

//  med.p  = regStep(inlet.p - outlet.p, inlet.p, outlet.p, dpTol) ;
  med.p  = sqrt(inlet.p * outlet.p) ;
  med.h  = inlet.h_outflow  ;
  med.Xi = inlet.Xi_outflow ;
   

end partialValve;


/*
// History
// Pressure difference drives the flow
//     Pressure used for state is regularized .. 
//     avoid sudden change during reversal 

//  None of these worked well with multiple valves connections
//  state = Medium.setState_phX( regStep( (inlet.p - outlet.p), inlet.p, outlet.p, 10),
//    			       inStream(outlet.h_outflow),
//    			       inStream(outlet.Xi_outflow));

// Smooth out state during flow reversal 
// This is more logical specification, and should work well.
// Now this specification is  working (1.9.3+dev (r25920))
// All valves in parallel with same source and sink converged (mmmmmmmmmmmmm)
// Ok my Adsorber is also happy
//  state = Medium.setSmoothState( (inlet.p - outlet.p), 
//                         Medium.setState_phX(inlet.p, inlet.h_outflow, inlet.Xi_outflow),
//                         Medium.setState_phX(outlet.p, outlet.h_outflow, outlet.Xi_outflow),
//    			10);  
// Crude smoothing
//  state = Medium.setState_phX((inlet.p+outlet.p)*0.5, 
//                              (inlet.h_outflow + outlet.h_outflow)*0.5,
//                              (outlet.Xi_outflow + inlet.Xi_outflow)*0.5);
// try me temporary
//  With this logic you get clock partitioning problem wiht multivalve setup
//  commenting it out. Using conditional expression on state fails...with translation
//  errors.   But works //  on indivdual arguments that state is a function of. I think
//  it is because state is a record rather than a simple variable. Need to raise this as
//  a bug. 
//     state =  if (inlet.p > outlet.p) then
//                    Medium.setState_phX(pp, outlet.h_outflow, outlet.Xi_outflow) 
//            else
//                  Medium.setState_phX(outlet.p, outlet.h_outflow, outlet.Xi_outflow) ;
//    end if;
   

// Lesson: avoid calling setState_phX .. then you are OK ... since it avoids
//           not so flash internal SingleNonLinearEquation Solver 
//    but then you need these three equations   
//        WELL ... this stopped working since updating the version
//             and "the Crude Smoothing works" I wonder why
//          Have they tinkered with something in Medium package ... since I reported
//     this problem ?


//        if (Medium.nS <> 1) then
//          state.X = (inlet.Xi_outflow + outlet.Xi_outflow)*0.5 ;
//        end if ;
//       state.p = (inlet.p+outlet.p)*0.5 ; 
//       Medium.specificEnthalpy(state) =  (outlet.h_outflow + inlet.h_outflow)*0.5;

//  state =  Medium.setState_phX(if (inlet.p > outlet.p) then inlet.p else outlet.p ,
//                               if (inlet.p > outlet.p) then inlet.h_outflow else outlet.h_outflow,
//                               if (inlet.p > outlet.p) then inlet.Xi_outflow else outlet.Xi_outflow) ;
//
//  Medium.ThermodynamicState state(p(start=1e5, fixed=false), 
//                                  T(start=300, fixed=false), 
//                                  X(start=Medium.reference_X)) ;
//  Medium.ThermodynamicState state1(p(start=1e5, fixed=false), 
 //                                 T(start=300, fixed=false), 
//                                  X(start=Medium.reference_X)) ;
//  Medium.ThermodynamicState state2(p(start=1e5, fixed=false), 
//                                  T(start=300, fixed=false), 
 //                                 X(start=Medium.reference_X)) ;
//  state1 = Medium.setState_phX(inlet.p,  inlet.h_outflow, inlet.Xi_outflow); 
 // state2 = Medium.setState_phX(outlet.p, outlet.h_outflow, outlet.Xi_outflow);
  //state = Medium.setSmoothState( 2*(inlet.p - outlet.p)/(inlet.p + outlet.p), state1, state2, 0.01);  
  // state = Medium.setSmoothState( (inlet.p - outlet.p), state1, state2, 100);  
*/
/*
    state.p = inlet.p ;
    state.X = if Medium.reducedX then
                cat(1, inlet.Xi_outflow, {1 - sum(inlet.Xi_outflow)}) 
              else
                inlet.Xi_outflow ;
     inlet.h_outflow = Medium.specificEnthalpy(state);  // indirect way of making state.T as tear variable
                                                      //  without ever invoking T_hX 
//  inlet.h_outflow = h_TX(state.T, state.X) ;        // this won't work .... Ptex error 
*/

/*
  state = Medium.setState_phX(inlet.p,  inlet.h_outflow, inlet.Xi_outflow); 
  state = Medium.setSmoothState( (inlet.p - outlet.p), 
                        Medium.setState_phX(inlet.p, inlet.h_outflow, inlet.Xi_outflow),
                       Medium.setState_phX(outlet.p, outlet.h_outflow, outlet.Xi_outflow),
 			10);  
   med.p = inlet.p ; med.T=350 ; med.X=if (Medium.reducedX) then
                                           cat(1, inlet.Xi_outflow, {1 - sum(inlet.Xi_outflow)})
                                       else
                                          inlet.Xi_outflow ;




   state1 = Medium.setState_phX(inlet.p,  inlet.h_outflow, inlet.Xi_outflow); 
   state2 = Medium.setState_phX(outlet.p, outlet.h_outflow, outlet.Xi_outflow);
   state = Medium.setSmoothState( (inlet.p - outlet.p) , state1, state2, dpTol);  

   state = Medium.setSmoothState((inlet.p - outlet.p), 
                                   Medium.setState_phX(inlet.p,  inlet.h_outflow, inlet.Xi_outflow), 
                                   Medium.setState_phX(outlet.p, outlet.h_outflow, outlet.Xi_outflow),
                                   dpTol);
*/

