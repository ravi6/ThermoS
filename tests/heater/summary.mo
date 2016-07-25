package Summary "Uop Summary Package"

import Modelica.Utilities.Strings.isEqual;
import OpenModelica.Scripting.* ;

function genReport

/* Generates report for a given 
   equipment and sample time from an
   existing results file
*/

     input String  resFile ;
     input String  uopName ;
     input Real    sTime   ;
    
     algorithm

/* Now populate report configurations for all types of Uops */
      //  if (isEqual(uopType, "heater", false)) then
     //print(stringTypeName(uopName)); 

         print("=============<< Heater Summary >>=============\n"  +
                   "Tag: " + uopName + 
                   "\tSample: " + String(sTime) +
                   "\nInlet Conditions:"   +
                   "\nFlow rate (kg/s) ="  + String(val(stringVariableName(uopName+".inlet.m_flow"), sTime, resFile)) +
                   "\nExit Conditions:"    + 
                   "\nFluid Temp (C) ="    + String(val(stringVariableName(uopName+".Tf"), sTime, resFile)) +
                   "\tWall Temp (C) ="     + String(val(stringVariableName(uopName+".Tw"), sTime, resFile)) +
                   "\nQ_ew (kW) = "        + String(val(stringVariableName(uopName+".Q_ew"), sTime, resFile)*1e-3) +
                   "\tQ_wf (kW) = "        + String(val(stringVariableName(uopName+".Q_wf"), sTime, resFile)*1e-3) +
                   "\tCp (J/kgC) = "       + String(val(stringVariableName(uopName+".Cp"), sTime, resFile)) +
                   "\n==========================================" 
               );
       // end if;
end genReport;

end Summary;
