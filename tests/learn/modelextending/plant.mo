model plant  "A plant"
mylib.bus abus(speed=5, k=3);

initial equation
 abus.x = 0 ;
 abus.y = 0 ;
end plant;
