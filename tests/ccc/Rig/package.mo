//Author: Ravi Saripalli
package Rig   "All Models needed for Rig"
   import Modelica.Media.Interfaces.PartialMixtureMedium;
   import Modelica.Fluid.Interfaces.FluidPort;
   import Modelica.Fluid.Utilities.*;  // Some regularization functions

   import Modelica.Units.SI.* ;        // Embrace SI units
   type Percent = Real(unit="p", min=0 , max=100) ;
   type Fraction = Real(min=0, max=1.0) ;
end Rig;
