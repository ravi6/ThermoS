function gen_getVar 
/* Generate script lines that can
    retrieve a Variable value at a specific sample 
    from an existing results file
*/
    input String resFile ;
    input String uopName ;
    input String varName ;
    input Integer count  ;
    
    output Boolean      ans;

     algorithm
        ans:=OpenModelica.Scripting.writeFile("/tmp/rep.mos",
           "data" + " := readSimulationResult( " +
           "\"" + resFile + "\"," + 
           uopName + "." + varName + ", size=0);" + "\n" +
           varName + ":=" + "data[1," + String(count) + "]; \n " 
         , true);
end gen_getVar ;

function genReport
/* Generates report for a given 
   equipment and sample number from an
   existing results file
*/
import Modelica.Utilities.Strings.isEqual;

     input String resFile ;
     input String  uopName ;
     input String  uopType ;
     input Integer count   ;

     output Boolean ans ;
     output String  msg ;
    
     algorithm
      
        ans := OpenModelica.Scripting.writeFile("/tmp/rep.mos","",false); // clean file
        ans := OpenModelica.Scripting.writeFile("/tmp/rep.mos",
                           "name" + ":=" +  "\"" + uopName  + "\"" + "; \n " + 
                           "count" + ":=" +  String(count) + "; \n " 
                        , true); 

/* Now populate report configurations for all types of Uops */
        if (isEqual(uopType, "heater", false)) then
            ans := gen_getVar(resFile, uopName, "Tf", count) ;
            ans := gen_getVar(resFile, uopName, "Tw", count) ;

            ans := OpenModelica.Scripting.writeFile("/tmp/rep.mos",
                        "print( \"<<Heater Summary>>\" + \"\\n\" + 
                                \"Tag: \" + name + \"\\t\"    +
                                \"Sample: \" + String(count) + \"\\n\"    +
                                \"Fluid Temp =\"  + String(Tf)  +  \"\\n\" + 
                               \"Wall Temp = \" + String(Tw)  + \"\\n\"
                               );"
                        , true); 
        end if;
        
        msg := OpenModelica.Scripting.runScript("/tmp/rep.mos");
end genReport;
