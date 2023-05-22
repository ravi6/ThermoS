model plant
  import ThermoS.Types.*;
  import ThermoS.Media.MyGas;
  import ThermoS.Uops.Feed;
  import ThermoS.Uops.Valves.Valve;
  import ThermoS.Uops.Tanks.GasTank;
  import ThermoS.Uops.Reservoir;
  import ThermoS.Uops.Controller;
  constant Real[MyGas.nXi] Air = {0.79, 0.21};
  Feed supply(redeclare package Medium = MyGas);
  Valve valve(redeclare package Medium = MyGas, cv = (1000e-3/60)/sqrt(4e5));
  GasTank tank(redeclare package Medium = MyGas, vol = 0.2, Q_in = 0);
  Reservoir atm(redeclare package Medium = MyGas, p = 1e5, T = 300, Xi = Air);
  Controller pid(Kc = 1e4, Ti = 100, Td = 0, reverseActing = true, pvMin = 1e5, pvMax = 11e5, mvMin = 0, mvMax = 6000e-3/60);
initial algorithm
  tank.T := 300;
  tank.p := 1e5;
  tank.Xi := Air;
equation
  connect(supply.outlet, tank.inlet);
  connect(tank.outlet, valve.inlet);
  connect(valve.outlet, atm.port);
  pid.sp = 1e5 + 6e5*(1 - exp(-time/10));
  pid.pv = tank.p;
  pid.mv = tank.inlet.m_flow;
  supply.T = 300;
  supply.Xi = fill(1.0/MyGas.nS, MyGas.nXi);
  valve.po = 50;
end plant;

package ThermoS
  package Uops
    package Valves
      partial model partialValve
        import Modelica.Media.IdealGases.Common.MixtureGasNasa.h_TX;
        replaceable package Medium = PartialMixtureMedium;
        FluidPort inlet(redeclare package Medium = Medium);
        FluidPort outlet(redeclare package Medium = Medium);
        parameter Real cv = 1.0/sqrt(1e5);
        parameter Real dpTol = 100;
        Medium.BaseProperties med;
      equation
        outlet.h_outflow = inStream(inlet.h_outflow);
        inlet.h_outflow = inStream(outlet.h_outflow);
        inlet.Xi_outflow = inStream(outlet.Xi_outflow);
        outlet.Xi_outflow = inStream(inlet.Xi_outflow);
        inlet.m_flow + outlet.m_flow = 0;
        med.p = max(inlet.p, outlet.p);
        med.h = inlet.h_outflow;
        med.Xi = inlet.Xi_outflow;
      end partialValve;

      model Valve
        extends ThermoS.Uops.Valves.partialValve;
        parameter Vchar vchar = Vchar.Linear;
        parameter Fraction pratChoke = 0.5;
        parameter Boolean Compressible = true;
        Percent po(start = 0.0);
        Fraction charF(start = 1.0);
        Fraction prat(start = 1.0);
      equation
        if (vchar == Vchar.Linear) then
          charF = po/100;
        elseif (vchar == Vchar.FastActing) then
          charF = (po/100)^0.5;
        elseif (vchar == Vchar.EquiPercent) then
          charF = 35^(po/100 - 1);
        end if;
        prat = min(inlet.p, outlet.p)/max(inlet.p, outlet.p);
        if (Compressible) then
          inlet.m_flow = noEvent(if (prat > 1 or prat < 1) then cv*max(0, charF)*sqrt(max(0, med.d))*sqrt(max(inlet.p, outlet.p))*sign(inlet.p - outlet.p)*regRoot(1 - max(prat, 0.5), dpTol) else 0);
        else
          inlet.m_flow = cv*charF*sqrt(med.d*inlet.p)*regRoot2(1 - outlet.p/inlet.p, dpTol);
        end if;
      end Valve;

      type Vchar = enumeration(Linear, FastActing, EquiPercent);
    end Valves;

    package Tanks
      model GasTank
        replaceable package Medium = PartialMixtureMedium;
        FluidPort inlet(redeclare package Medium = Medium);
        FluidPort outlet(redeclare package Medium = Medium);
        parameter Volume vol = 10;
        parameter EnthalpyFlowRate Q_in = 0;
        Mass m;
        Medium.Temperature T;
        Medium.AbsolutePressure p;
        Medium.MassFraction[Medium.nXi] Xi;
        Medium.SpecificEnthalpy h;
        Medium.BaseProperties medium;
      equation
        medium.T = T;
        medium.p = p;
        medium.Xi = Xi;
        medium.h = h;
        m = medium.d*vol;
        der(m) = inlet.m_flow + outlet.m_flow;
        der(Xi*m) = actualStream(inlet.Xi_outflow)*inlet.m_flow + actualStream(outlet.Xi_outflow)*outlet.m_flow;
        der(m*h) = Q_in + inlet.m_flow*actualStream(inlet.h_outflow) + outlet.m_flow*actualStream(outlet.h_outflow) + vol*der(p);
        inlet.Xi_outflow = Xi;
        inlet.h_outflow = h;
        outlet.Xi_outflow = Xi;
        outlet.h_outflow = h;
        inlet.p = outlet.p;
        inlet.p = p;
      end GasTank;
    end Tanks;

    model Reservoir
      replaceable package Medium = PartialMixtureMedium;
      FluidPort port(redeclare package Medium = Medium);
      parameter Medium.AbsolutePressure p;
      parameter Medium.Temperature T;
      parameter Medium.MassFraction[Medium.nXi] Xi;
      Medium.BaseProperties medium;
    equation
      medium.T = T;
      medium.p = p;
      medium.Xi = Xi;
      port.h_outflow = medium.h;
      port.Xi_outflow = Xi;
      port.p = p;
    end Reservoir;

    package Interfaces end Interfaces;

    model Feed
      replaceable package Medium = PartialMixtureMedium;
      FluidPort outlet(redeclare package Medium = Medium);
      Medium.MassFlowRate mdot(min = 0);
      Medium.Temperature T;
      Medium.MassFraction[Medium.nXi] Xi;
      Medium.BaseProperties medium;
    equation
      medium.p = outlet.p;
      medium.T = T;
      medium.Xi = Xi;
      outlet.h_outflow = medium.h;
      outlet.m_flow = -mdot;
      outlet.Xi_outflow = Xi;
    end Feed;

    model Controller
      Real Kc;
      Real Ti;
      Real Td;
      parameter Boolean reverseActing = true;
      parameter Real mvMin = 0;
      parameter Real mvMax = 100;
      parameter Real pvMin = 0;
      parameter Real pvMax = 100;
      parameter Integer action = if reverseActing then 1 else -1;
      Real sp(min = pvMin, max = pvMax);
      Real pv;
      Real mv(min = mvMin, max = mvMax);
      Real op(min = mvMin, max = mvMax);
      Real err;
      Real intErr;
    initial equation
      intErr = 0;
    equation
      err = (sp - pv);
      der(intErr) = noEvent(if mv < mvMin and err < 0 or mv > mvMax and err > 0 then 0 else err);
      op = action*((mvMax - mvMin)/(pvMax - pvMin))*Kc*(err + intErr/Ti + Td*der(err)) + mvMin;
      mv = noEvent(if op < mvMin then mvMin elseif op > mvMax then mvMax else op);
    end Controller;

    import Modelica.Media.Interfaces.PartialMixtureMedium;
    import Modelica.Fluid.Interfaces.FluidPort;
    import Modelica.Fluid.Utilities.*;
    import Modelica.Units.SI.*;
    import ThermoS.Types.*;
  end Uops;

  package Types
    type Percent = Real(unit = "p", min = 0, max = 100);
    type Fraction = Real(min = 0, max = 1.0);
  end Types;

  package Media
    package MyGas
      import Modelica.Media.IdealGases.Common.MixtureGasNasa;
      import Modelica.Media.IdealGases.Common.SingleGasesData;
      import Modelica.Media.IdealGases.Common.FluidData;
      extends MixtureGasNasa(data = {SingleGasesData.N2, // note data is of type  DataRecord[:]
      SingleGasesData.O2, SingleGasesData.CO2}, fluidConstants = {FluidData.N2, FluidData.O2, FluidData.CO2}, substanceNames = {"Nitrogen", "Oxygen", "CarbonDioxide"}, reducedX = true, reference_X = {0.7, 0.2, 0.1}, Density(start = 1, nominal = 1), AbsolutePressure(start = 1e5, min = 1e3, max = 50e5, nominal = 1e5), Temperature(start = 300, min = 200, max = 2000, nominal = 300));
    end MyGas;
  end Media;
end ThermoS;

model plant_total
  extends plant;
end plant_total;
