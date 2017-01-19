model plant  
  Real ymin;
  Real ymax;
  Real x(start = -1);
initial equation
  x = -10;
equation
  der(x) = 1;
  ymin = ThermoS.Math.sMin(x, 5, 0.01);
  ymax = ThermoS.Math.sMax(x, 3, 0.01);
end plant;

package ThermoS  "A Modelica Package for Process Simulations" 
  package Math  "All of Math Functionality here " 
    function sMin  
      input Real x;
      input Real y;
      input Real eps;
      output Real z;
    algorithm
      z := -(sqrt((x - y) ^ 2 + eps ^ 2) - (x - y)) / 2;
      annotation(Inline = "true", smoothOrder = 2); 
    end sMin;

    function sMax  
      input Real x;
      input Real y;
      input Real eps;
      output Real z;
    algorithm
      z := (sqrt((x - y) ^ 2 + eps ^ 2) + x - y) / 2;
      annotation(Inline = "true", smoothOrder = 2); 
    end sMax;
  end Math;
end ThermoS;

model plant_total
  extends plant;
end plant_total;
