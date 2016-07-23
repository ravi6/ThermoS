//Author: Ravi Saripalli
within ThermoS.Uops;
package Valves   "Contains Valves and its interfaces "  

        type Vchar = enumeration ( Linear         "Linear Valve",
                                   FastActing     "Fast Acting Valve",
                                   EquiPercent    "Equi-Percent Valve"
                                 ) "Enumeration Defining Valve Behaviour" ; 
end Valves;
