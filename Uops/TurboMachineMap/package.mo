within ThermoS.Uops;
/* Version 1.2
   Last Modified 13th March 2013
   Author: Ravi Saripalli
 * History  
    * Added Choking and Surge Limiting Functions
    * Added more functions to permit Turbine calcs
*/

package TurboMachineMap "Implements Compressor/Turbine Map functionality"
import Modelica.Blocks.Types.ExternalCombiTable2D;
import Modelica.Blocks.Types.ExternalCombiTable1D;
import Modelica.Blocks.Types.Smoothness;
/* Mc as function of beta at different speeds */
constant Real[:, :] mcTable=[0.0, 1, 2, 3, 4, 5, 6; 
			     0.0, 2, 3, 4, 5, 6, 7;
 			     0.5, 3, 4, 5, 6, 7, 8; 
			     1.0, 4, 5, 6, 7, 8, 9];
/* Prat as function of beta at different speeds */
constant Real[:, :] prTable=[0, 1, 2, 3, 4, 5, 6;
			     0.0, 3, 4, 5, 6, 7, 8; 
                             0.5, 2, 3, 4, 5, 6, 7; 
			     1.0, 1, 2, 3, 4, 5, 6];
/* Eff as function of beta at different speeds */
constant Real[:, :] effTable=[0, 1, 2, 3, 4, 5, 6;
                              1, 1, 1, 1, 1, 1, 1];
constant Temperature 		Tref = 300; // Kelvin
constant AbsolutePressure 	Pref = 1e5; // Pascals
record MapPoint // Holds Map Data at a specific operating condition
  Real 			beta ; 		// parameteric line value
  Real 			Nc ; 		// Corrected Speed
  MassFlowRate 		Mc ; 		// Corrected Mass Flow Rate
  Real 			pratio ;	 // Pressure Ratio
  Real 			Eff ; 		// Adiabatic Efficiency
  end 			MapPoint ;
/****************/
function Intpol2D
/****************/
  input Real x;
  input Real y;
  input ExternalCombiTable2D TableID;
  output Real z;
/*  Now interpolate  */
external"C" z = ModelicaStandardTables_CombiTable2D_getValue(
      TableID, x, y) annotation (Library={"ModelicaStandardTables"});
  end Intpol2D;
/****************/
function Intpol1D
/****************/
  input Real x;
  input ExternalCombiTable1D TableID;
  output Real y;
/*  Now interpolate  nb. We interpolate on only one column*/
external"C" y = ModelicaStandardTables_CombiTable1D_getValue( TableID, 1, x) 
                annotation (Library={"ModelicaStandardTables"});
  end Intpol1D;
/**************/
function getMc
/**************/
/* Gets Corrected Mass flow from the map */
  input Real beta;
  input Real Nc;
  output Real mc;
protected
  ExternalCombiTable2D TableID=ExternalCombiTable2D(
        				"NoName", "NoName", mcTable,
        				Smoothness.LinearSegments);
algorithm
     mc := Intpol2D( beta, Nc, TableID);
end getMc;
/*****************/
function getPratio
/*****************/
/* Gets Pressure Ratio from the map */
  input Real beta;
  input Real Nc;
  output Real pr;
protected
  ExternalCombiTable2D TableID=ExternalCombiTable2D( 
			       "NoName", "NoName", prTable,
        			Smoothness.LinearSegments);
algorithm
  pr := Intpol2D( beta, Nc, TableID);
end getPratio;
/*****************/
function getEff
/*****************/
/* Gets Adiabatic Efficiency from the map */
  input Real beta;
  input Real Nc;
  output Real eff;
protected
  ExternalCombiTable2D TableID=ExternalCombiTable2D(
        			"NoName", "NoName", effTable,
        			Smoothness.LinearSegments);
algorithm
  eff := Intpol2D( beta, Nc, TableID);
end getEff;
/************************/
function fnChokedFlow
/***********************/
/* Choked Mass Flow Rate */
  input AngularVelocity N;
  input Temperature T_in;
  input AbsolutePressure P_in;
  output MassFlowRate mflow_choke;
protected Real MCorFactor;
protected Real NCorFactor;
protected Real Nc;
  algorithm
  NCorFactor := 1.0/sqrt(T_in/Tref);
  MCorFactor := 1.0/(NCorFactor*(P_in/Pref));
  Nc := N*NCorFactor;
  mflow_choke := getMc(1.0, Nc)/MCorFactor;
end fnChokedFlow;
/************************/
function fnMc
/***********************/
 //Determine Corrrect Mass Flow
  input Temperature T_in;
  input AbsolutePressure P_in;
  input MassFlowRate mflow_in;
  output MassFlowRate   Mc;
protected Real MCorFactor;
protected Real NCorFactor;
  algorithm
  NCorFactor := 1.0/sqrt(T_in/Tref);
  MCorFactor := 1.0/(NCorFactor*(P_in/Pref));
  Mc := mflow_in*MCorFactor;
end fnMc; 
/************************/
function fnNc
/***********************/
 //Determine Corrrected Speed
  input Temperature T_in;
  input AbsolutePressure P_in;
  input AngularVelocity N;
  output Real 		Nc;
protected Real MCorFactor;
protected Real NCorFactor;
  algorithm
  NCorFactor := 1.0/sqrt(T_in/Tref);
  MCorFactor := 1.0/(NCorFactor*(P_in/Pref));
  Nc := N*NCorFactor;
end fnNc; 
end TurboMachineMap;
