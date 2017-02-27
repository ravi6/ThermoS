model fo  
  Real y;
  Real x;
  parameter Real tau = 10;
initial equation
  y = x;
equation
  tau * der(y) + y = x;
end fo;

model plant  
  fo filter;
equation
  filter.x = sin(time);
end plant;

model plant_total
  extends plant;
end plant_total;
