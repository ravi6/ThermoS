within ThermoS.Media;

package JP8 "Aviation Fuel JP8"

extends Modelica.Media.Incompressible.TableBased (
                 mediumName = "Aviation Fuel JP 8",
                 T_min = Modelica.SIunits.Conversions.NonSIunits.from_degC(-50), 
                 T_max = Modelica.SIunits.Conversions.NonSIunits.from_degC(230), 
                 TinK = false,      // Table Temperatures are in C
                 T0 = 273.15, 
                  reducedX = true,
              tableDensity = [-40.0, 850; 0.0, 822 ; 20.0, 808 ; 90.0, 760 ],
              tableHeatCapacity = [31.3, 2000; 60, 2125; 100, 2300; 183, 2660],
              tableConductivity = [-6.8, 0.12; 49.3, 0.101; 118, 0.098; 219.5, 0.0803],
              tableViscosity = [20, 1.368e-03; 30, 1.156e-03; 40, 9.933e-04;
                                50, 8.657e-04; 60, 7.632e-04; 70, 6.785e-04;
                                80, 6.050e-04; 90, 5.451e-04; 100, 4.939e-04],
              tableVaporPressure = [-40,  5.81; -20,  27.15; 0,  101.3; 20,  315.7;
                                     40,  851 ;  60,  2037 ; 80, 4415; 100,  8811;
                                     120, 16387; 140, 28701; 160, 47733; 180, 75900; 200, 116048]);

/* Property tables derived from 
(1) Thermodynamic, Transport, and Chemical Properties of Reference JP-8
Thomas J. Bruno et. al, NISTIR 6659,Physical and Chemical Properties Division
National Institute of Standards and Technology, 325 Broadway, Boulder, CO 80305-3337
July 2010 
(2) Aviation Fuel Properties, CRC Report No. 530
Coordinating Research Council, Inc. Atlanta, Georgia
Handbook of Aviation Fuel Properties, prepared by CRC Inc.
*/
end JP8;
