package Summary  "UOP Summary package"
import OpenModelica.Scripting.* ;
import Modelica.Utilities.Strings.isEqual;


//constant String CR = "\n" ; 

function qs
// Generate a quoted string
   input String strin   ;
   output String strout ;

algorithm
   strout := "\"" + strin + "\"" ;
end qs ;



function gen_getVar 
/* Generate script lines that can
    retrieve a Variable value at a specific sample Time 
    from an existing results file and assign it to
    a new variable name (Alias ... this avoids dot problem)
*/
    input String resFile ;
    input String uopName ;
    input String varName ;
    input String varAlias ;
    input Real   sTime ;        // Sample time in seconds
    output Boolean      ans;

    protected String str ;
    

     algorithm
        str := varAlias + ":=" + 
               "val(" + 
               uopName + "." + varName + "," + String(sTime) +
               "," +  qs(resFile) + "); \n"     ;
        ans := writeFile("/tmp/rep.mos", str, true);
        
end gen_getVar ;

function genReport
/* Generates report for a given 
   equipment and sample number from an
   existing results file
*/

     input String resFile ;
     input String  uopName ;
     input String  uopType ;
     input Real   sTime ;        // Sample time in seconds

     output Boolean ans ;
     output String  msg ;
    
     algorithm
      
        ans := writeFile("/tmp/rep.mos","",false); // clean file
        ans := writeFile("/tmp/rep.mos",
                         "name" + ":=" +  "\"" +
                         uopName  + "\"" + "; \n " + 
                         "sTime" + ":=" +  String(sTime) + "; \n " 
                         , true); 

/* Now populate report configurations for all types of Uops */
        if (isEqual(uopType, "heater", false)) then
            ans := gen_getVar(resFile, uopName, "Tf", "Tf", sTime) ;
            ans := gen_getVar(resFile, uopName, "Tw", "Tw", sTime) ;
            ans := gen_getVar(resFile, uopName, "Q_ew", "Q_ew", sTime) ;
            ans := gen_getVar(resFile, uopName, "Q_wf",  "Q_wf",sTime) ;
            ans := gen_getVar(resFile, uopName, "inlet.m_flow", "mflow", sTime) ;

            ans := OpenModelica.Scripting.writeFile("/tmp/rep.mos",
                        "print( \"<<Heater Summary>>\" + \"\\n\" + 
                                \"Tag: \" + name + \"\\t\"    +
                                \"Sample Time: \" + String(sTime) + \"\\n\"    +
                                \"Inlet Conditions:\" +    \"\\n\" + 
                                \"flow rate (kg/s) =\"  + String(mflow)  +  \"\\n\" + 
                                \"Exit Conditions:\" +    \"\\n\" + 
                                \"Fluid Temp (C) =\"  + String(Tf)  +  \"\\t\" + 
                                \"Wall Temp (C) =\"  + String(Tw)  +  \"\\n\" + 
                                \"Q_ew (kW) = \" + String(Q_ew/1e3)  + \"\\t\" +
                                \"Q_wf (kW) = \" + String(Q_ew/1e3)  + \"\\n\"
                               );"
                        , true); 
        end if;
        
        msg := OpenModelica.Scripting.runScript("/tmp/rep.mos");
end genReport;

/* Looks very complicated for no reason ...
          but the scripting engine does not
   permit dynamic variable loading through the use of 
readSimulationResults calls... the only recourse is to run it as
a separate script ...
*/
end Summary ;
