package Summary "Uop Summary Package"

import Modelica.Utilities.Strings.isEqual;
import OpenModelica.Scripting.* ;

function genReport

/* Generates report for a given 
   equipment and sample time from an
   existing results file
*/
     input String  uopName ;
     input String  resFile ;
     input Real    sTime   ;

    protected Real data[:, :];
    algorithm 
         readSimulationResult("plant_res.csv", stringVariableName("time"));
         print (getErrorString());
         print (String(size(data,1))) ;
         print ("  " + String(size(data,2))) ;
       
/*
      print("=============<< Heater Summary >>=============");
      print("Tag: " + uopName + "\tSample: " + String(sTime));
      print("\nInlet Conditions:");
      print("Flow rate (kg/s) = " + String(readVar(uopName+".inlet.m_flow", sTime, resFile)));
      print("\nExit Conditions:");
      print("Fluid Temp (C) = " + String(readVar(uopName+".Tf", sTime, resFile)) +
	    "\tWall Temp (C) = " + String(readVar(uopName+".Tw", sTime, resFile)));
      print("Q_ew (kW) = " + String(readVar(uopName+".Q_ew", sTime, resFile)*1e-3) +
	    "\tQ_wf (kW) = " + String(readVar(uopName+".Q_wf", sTime, resFile)*1e-3) +
	    "\tCp (J/kgC) = " + String(readVar(uopName+".Cp", sTime, resFile)));
      print("==========================================");
*/
end genReport;


function readVar 
  input String varName;
  input Real t;
  input String fileName;
  output Real var;

protected
  Real data[:, :];
  Integer n;

algorithm
  data := readSimulationResult(fileName, strVariable(varName));
  print(varName + "\n") ;
  n := size(data, 1);

  // Linear interpolation
  for i in 2:n loop
    if data[i,1] >= t then
      var := data[i-1,2] + (t - data[i-1,1]) *
                (data[i,2] - data[i-1,2]) /
                (data[i,1] - data[i-1,1]);
      return;
    end if;
  end for;

  // If time beyond last point, return last value
  var := data[n,2];
end readVar;

end Summary;
