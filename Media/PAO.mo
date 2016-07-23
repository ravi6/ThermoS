within ThermoS.Media;
package PAO "Poly Aliphatic Olefine Oil"
extends Modelica.Media.Incompressible.TableBased (
                 mediumName = "PAO Oil",
                 T_min = Modelica.SIunits.Conversions.NonSIunits.from_degC(-50), 
                 T_max = Modelica.SIunits.Conversions.NonSIunits.from_degC(200), 
                 TinK = false,      // Table Temperatures are in C
                 T0 = 273.15, 
              tableDensity = [4.7, 800;  41.1, 779;  71.1, 759;
                              95.6, 742;  123.5, 723;  155.7, 702;    194, 673],
              tableHeatCapacity = [-0.6, 2103; 42, 2280; 95, 2466; 148, 2648;
                                    204, 2835; 239, 2964],
              tableConductivity = [60.4, 0.1444; 101.7, 0.1387; 161.4, 0.1313;
                                   207, 0.1255;  242.6, 0.1212], 
              tableViscosity = [-49.2, 0.497; -46.8, 0.3563; -41.5, 0.2112; -26.9, 0.1264;
                                -1.36, 0.0436; 19.3, 0.0106; 40.6, 0.0079; 79.4, 0.00482;
                                100 ,0.002164],
              tableVaporPressure = [0, 500; 20, 1.9e3; 40, 5.3e3; 
                                    60, 16e3; 80, 37e3; 100, 80e3]);  // wrong input ...fix it
/* Property tables derived from 
JOURNAL OF SYNTHETIC LUBRICATION
J. Synthetic Lubrication 2007; 24: 7790
Published online 2 March 2007 in Wiley InterScience
(www.interscience.wiley.com) DOI: 10.1002/jsl.30
Advanced cooling technology for rotors of high-power
low-duty cycle generators using polyalphaolefins
Ahmad K. Sleiti
*/
end PAO;
