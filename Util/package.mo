within ThermoS;
package Util "General Utility functions"
/***********************************************************/
/* Useful String Manipulation Functions */
/* These exist in Modelica.Utilities.Strings .. 
      but just an illustration of using external calls */
// Added DataSave sub_package to permit clean MatLab readable 
//   files that can have variables of all types
/***********************************************************/
import OpenModelica.Scripting.* ;
import Modelica.Utilities.Strings.isEqual;

/*==================*/
 function strlen 
/*==================*/
    input String str; 
    output Integer len; 
    external "C" len=ModelicaStrings_length(str);
 end strlen;
/*==================*/
function substr 
/*==================*/
    input String str; 
    input Integer i1,i2; 
    output String out; 
    external "C" out=ModelicaStrings_substring(str,i1,i2);
 end substr;
/* Other functions I need for pretty printing 
   Error Checking in Simulation etc.
/*==================*/
function decimalS
/*==================*/
   input Real    x ;
   input Integer n ;
   output String   s ;
 algorithm
    s := String(floor(x*(10^n)) /  (10^n)) ;      
end decimalS;

/*=============================*/
function strVec
/*=============================*/
// Stringify an Array
  input  Real     v[:]  ;
  output String   s     ;

 algorithm
    s := "[" ;
    for i in 1:size(v,1)-1 loop
      s := s + String(v[i]) + ", "; 
    end for;
      s := s + String(v[size(v,1)]) + "] "; 
end strVec;      


/*==================*/
function ldFile
/*==================*/
/* Load file and bail out if unsuccessful */
	input String fname;
	output Boolean result;
 algorithm
        result := loadFile(fname) ;
	if (result) then
	       print(fname + " is Loaded\n");
        else
	       print(fname + " failed to load\n");
	       exit(1) ;
	end if;
end ldFile;


/*=================*/
function getInfo
/*================*/
  // Some useful info on Modelica Options
/* See what kind of TearingMethods Exist */

protected 
   String methods[:] ;
   String src[:] ;

algorithm

    (methods, src):=getAvailableTearingMethods(); 
    print("\n** Available Tearing Methods\n");
    for k in (1:size(methods,1)) loop
        print(methods[k] + "\t\t::\t " + src[k] + "\n") ;
    end for;

    print("\n** Available IndexReduction Methods\n");
    (methods, src):=getAvailableIndexReductionMethods(); 
    for k in (1:size(methods,1)) loop
        print(methods[k] + "\t\t::\t" + src[k] + "\n") ;
    end for;

    (methods, src):=getAvailableMatchingAlgorithms(); 
    print("\n** Available Matching Algorithms\n");
    for k in (1:size(methods,1)) loop
        print(methods[k] + "\t\t::\t" + src[k] + "\n") ;
    end for;


    print("\n ** Solver Settings  ** \n");
    print("TearingMethod          : " + getTearingMethod() + "\n");
    print("Index Reduction Method : " + getIndexReductionMethod() + "\n");
    print("Matchin Alogrithm      : " + getMatchingAlgorithm() + "\n");

end getInfo;



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

function printVars
// print a bunch of variables quickly

     input String  resFile ;
     input String[:]   varNames ;
     input Real    sTime   ;
     algorithm
                for k in 1:size(varNames,1) loop
                   print( varNames[k] +  " = " 
                    + String(val(stringVariableName(varNames[k]), sTime, resFile)) 
                    + "\n" );
                 end for;
end printVars;

package DataSave  "Contains Matlab Friendly Save calls"

    replaceable constant  String outFile = "data" ;
    replaceable constant  String resFile = "plant_res.mat" ;
    replaceable constant  Real   sTime = 0.0 ;

function saveVars
// print a bunch of variables quickly

     input String[:]   varNames ;

    protected 
          Boolean ans ;

     algorithm
            for k in 1:size(varNames,1) loop
                ans := writeFile(outFile, "data." + varNames[k] +  " = " 
                                  + String(val(stringVariableName(varNames[k]), sTime, resFile)) 
                                  + ";\n" , true);
            end for;
end saveVars;

function saveArray
// print a bunch of variables quickly

     input String   varName ;
     input Integer  M ;
     input Integer  N ;    // set this to zero for one Array


    protected 
      Boolean ans ;
      String str ;

     algorithm
             if ( N <> 0 ) then
                ans := writeFile(outFile, "data." + varName + "=" + "[" + "\n",true);
                for i in 1:M loop
                  for j in 1:N-1 loop
                    str := varName + "[" +  String(i) + "," + String(j) + "]" ;
                    ans := writeFile(outFile, String(val(stringVariableName(str), sTime, resFile)),true) ;
                    ans := writeFile(outFile, ", ",true);
                  end for;
                    str := varName + "[" +  String(i) + "," + String(N) + "]" ;
                    ans := writeFile(outFile, String(val(stringVariableName(str), sTime, resFile)),true) ;
                    ans := writeFile(outFile, ";",true);
                end for;
                    ans := writeFile(outFile, "];\n" ,true);

             else
                ans := writeFile(outFile, "data." + varName + "=" + "[" + "\n",true);
                for i in 1:M-1 loop
                    str := varName + "[" +  String(i)  + "]" ;
                    ans := writeFile(outFile, String(val(stringVariableName(str), sTime, resFile)),true) ;
                    ans := writeFile(outFile, ", ",true);
                end for;
                    str := varName + "[" +  String(M) + "]" ;
                    ans := writeFile(outFile, String(val(stringVariableName(str), sTime, resFile)),true) ;
                    ans := writeFile(outFile, "];\n" ,true);
             end if;

end saveArray;

end DataSave;  



end Util;
