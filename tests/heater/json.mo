function jsonResults
/* Generate a json data set of results
   at a specified time
   from a OpenModelica result file

     Author: Ravi Saripalli
     Date  : 30th July 2014
*/
import OpenModelica.Scripting.* ;

    input  String   resFile ;
    input  Real     sTime ;        // Sample time in seconds
    output Boolean  ans;

    protected String        str ;
    protected String[:]     varNames ;
    protected Integer       nvars ;
    protected Real          value; 

     algorithm

        varNames := readSimulationResultVars(resFile) ; 
        nvars    := 10 ; //size(varNames,1) ;  /* size fn. does not returns 0 instead of correct number */
                             

      // Construct the JSON string of results
        str := "{" ;
        for  i in 1:nvars-1 loop
             value := val(varNames[i], sTime, resFile) ;           
             str := str +  "\"" + varNames[i] + "\": " + String(value) + ",\n" ; 
        end for; 
             value := val(varNames[nvars], sTime, resFile) ;
             str := str +  "\"" + varNames[nvars] + "\": " + String(value) + "}\n" ; 
            
             ans := writeFile("/tmp/data.json", "", false); // clean file
             ans := writeFile("/tmp/data.json", str, true); // dump the string
       
end jsonResults ;
