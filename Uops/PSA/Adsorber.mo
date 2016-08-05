within ThermoS.Uops.PSA;
model Adsorber extends partAdsorber ; 

/*
Author: Ravi Saripalli
Version:  1.0
Date   : 3rd Jul 2015
Last Modfied :
*/

/* A full blown PSA Adsorber with two ports 
   Note:  -  Although ports are named as inlet and outlet they are just notional
             the unit should handle flow in any direction gracefully
          -  No energy balance in this model . The bed temperature is arbitrarily
             set to a constant value. (iso thermal) 
          -  We assume that Medium uses structurally independant composition spec.
             that is only Nc-1 compositions are evaluted and And the composition of
             the last element is = 1.0 - sum(yi)
*/

  import Modelica.Fluid.Interfaces.* ;
  import Modelica.Media.Interfaces.* ;

  bedParamsRec   bedParams ;
  // A packedbed adsorber with two ports (we want Medium to use reducedX=false)
  replaceable package Medium = PartialMixtureMedium ;

  // Specify that our Medium is used in/outlet
  FluidPort inlet(redeclare package Medium = Medium, Xi_outflow(start=Medium.reference_X[1:Medium.nXi]));
  FluidPort outlet(redeclare package Medium = Medium, Xi_outflow(start=Medium.reference_X[1:Medium.nXi]));

  Medium.ThermodynamicState     inlet_inState, outlet_inState ;
  Medium.ThermodynamicState     inlet_outState, outlet_outState; 

  constant Medium.MolarMass Mw[Medium.nX]  = Medium.fluidConstants[:].molarMass * 1000;

equation

    sum(yin_in)  = 1.0 ;
    sum(yin_out) = 1.0 ;

// Boundary Conditions
// ********************************************************************
// Inlet Mole fractions of fluid when flow is into the device at the two ports
// Note: yin_in/out are of size 1:nX(nS)
// ********************************************************************
    for i in 1:Medium.nXi loop
     inStream(inlet.Xi_outflow[i]) * sum (yin_in .* Mw)   =   yin_in[i] * Mw[i] ;
     inStream(outlet.Xi_outflow[i]) * sum (yin_out .* Mw) =  yin_out[i] * Mw[i] ;
    end for;

// Pressure Boundary Conditions
  p_in =  inlet.p / bedParams.Pref  ;
  p_out = outlet.p / bedParams.Pref ;

// ********************************************************************
//  The following is to provided connectivity of adsorber to external world
//   with the two fluid ports 
// ********************************************************************

    for i in 1:Medium.nXi loop
       inlet.Xi_outflow[i] * sum (y[N-1, :] .* Mw) =  y[N-1, i] * Mw[i] ; 
      outlet.Xi_outflow[i] * sum (y[N, :]   .* Mw) =  y[N, i]   * Mw[i] ; 
    end for;

// Connect mass flows to velocity
  
 // Avoiding T_hX calls (assuming reducedX)
//==========================================================================
  inlet_inState.p = inlet.p;
  inlet_inState.X =  cat(1, inStream(inlet.Xi_outflow), {1 - sum(inStream(inlet.Xi_outflow))});
  inStream(inlet.h_outflow) = Medium.specificEnthalpy(inlet_inState);

  outlet_inState.p =  outlet.p;
  outlet_inState.X =  cat(1, inStream(outlet.Xi_outflow), {1 - sum(inStream(outlet.Xi_outflow))});
  inStream(outlet.h_outflow) = Medium.specificEnthalpy(outlet_inState);
//==========================================================================

// These will call T_hX ... but very terse compared to above
//  inlet_inState =   Medium.setState_phX(inlet.p, inStream(inlet.h_outflow), inStream(inlet.Xi_outflow))  ;
//  outlet_inState =  Medium.setState_phX(outlet.p, inStream(outlet.h_outflow), inStream(outlet.Xi_outflow))  ;

  inlet_outState  = Medium.setState_pTX(inlet.p, bedParams.Tbed, inlet.Xi_outflow);
  outlet_outState = Medium.setState_pTX(outlet.p, bedParams.Tbed, outlet.Xi_outflow);


  inlet.m_flow = sign(u[N-1]) * bedParams.Uref * bedParams.csArea 
                  * (max(u[N-1], 0) * Medium.density(inlet_inState) +  
                      max(-u[N-1], 0) * Medium.density(inlet_outState)) ;


// Note flow convetion dictates that outflows are negative 
  outlet.m_flow = - sign(u[N])  * bedParams.Uref * bedParams.csArea 
                  * (max(-u[N], 0) * Medium.density(outlet_inState)  +
                     max(u[N], 0) *  Medium.density(outlet_outState));

// Set Enthalpys of outflow streams (that would result if flow is out of the device)
  inlet.h_outflow = Medium.specificEnthalpy(inlet_outState) ;
  outlet.h_outflow = Medium.specificEnthalpy(outlet_outState) ;

end Adsorber;

 /* Avoiding T_hX calls
  inlet_inState.p = inlet.p;
  inlet_inState.X = if Medium.reducedX then  cat(1, inStream(inlet.Xi_outflow), {1 - sum(inStream(inlet.Xi_outflow))})
                    else inlet_inState.X ;
  inStream(inlet.h_outflow) = Medium.specificEnthalpy(inlet_inState);

  outlet_inState.p =  outlet.p;
  outlet_inState.X = if Medium.reducedX then  cat(1, inStream(outlet.Xi_outflow), {1 - sum(inStream(outlet.Xi_outflow))})
                    else outlet_inState.X ;
  inStream(outlet.h_outflow) = Medium.specificEnthalpy(outlet_inState);

//     inStream(inlet.Xi_outflow[i])  =   ( 1 / sum (yin_in .* Mw))    *   yin_in[i] * Mw[i] ;
//     inStream(outlet.Xi_outflow[i]) =   ( 1 / sum (yin_out .* Mw))  *   yin_out[i] * Mw[i] ;
 */

/*
  Temporarymove
  inlet.m_flow = u[N-1] * bedParams.Uref * bedParams.csArea 
                  * homotopy( actual = (regStep(u[N-1], Medium.density(inlet_inState), 
                                    Medium.density(inlet_outState), uTol)),
                              simplified = 1.0 ) ;

//                  * regStep(p_in-p[1], Medium.density(inlet_inState), Medium.density(inlet_outState)) ;

// Note flow convetion dictates that outflows are negative 
  outlet.m_flow = - u[N]  * bedParams.Uref * bedParams.csArea 
                  * homotopy( actual = (regStep(-u[N], Medium.density(outlet_inState), 
                                   Medium.density(outlet_outState), uTol)),
                              simplified = 1.0 ) ;

//                * regStep(p[N]-p_out, Medium.density(outlet_inState), Medium.density(outlet_outState)) ;
*/

//                  * regStep(p_in-p[1], Medium.density(inlet_inState), Medium.density(inlet_outState)) ;
//                * regStep(p[N]-p_out, Medium.density(outlet_inState), Medium.density(outlet_outState)) ;
