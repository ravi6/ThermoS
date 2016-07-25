model plant  "A plant"
//mylib.bus abus(speed=15);
//mylib.geep ageep(redeclare package v=mylib.p);
mylib.geep ageep(redeclare model v=mylib.vehicle(speed=3));

initial equation
// abus.x = 0 ;
 //abus.y = 0 ;
end plant;
