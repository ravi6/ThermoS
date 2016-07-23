within ThermoS.Uops;
model Pump
     extends PumpBasic ;

  /* Pump Model for incompressible fluids 

*/
/*
  replaceable package PCurves = ThermoS.Interfaces.PumpCurves ;

  AngularVelocity 	speed(min = 0.1)            ; //(avoid singularity near zero speed)

  equation
    // Shaft work that is transferred to fluid
        head = PCurves.get_head(Vdot, speed);
        eff  = PCurves.get_eff(Vdot, speed);
*/
end Pump;
