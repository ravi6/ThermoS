model plant
  import ThermoS.Types.*;
  import ThermoS.Media.MyGas;
  import ThermoS.Uops.Feed;
  import ThermoS.Uops.Valves.Valve;
  import ThermoS.Uops.Tanks.OpenTank;
  import ThermoS.Uops.Reservoir;
  constant Real[MyGas.nXi] Air = {0.79, 0.21};
  Feed supply(redeclare package Medium = MyGas);
  Valve valve(redeclare package Medium = MyGas);
  OpenTank tank(redeclare package Medium = MyGas, in_pos = 0.44);
  Reservoir lake(redeclare package Medium = MyGas, p = 1e5, T = 300, Xi = Air);
initial algorithm
  tank.pFull := 50;
  tank.Tf := 400;
equation
  supply.mdot = 5 + 4.95*sin(6*time);
  supply.T = 300;
  supply.Xi = fill(1.0/MyGas.nS, MyGas.nXi);
  connect(supply.outlet, tank.inlet);
  connect(tank.outlet, valve.inlet);
  connect(valve.outlet, lake.port);
  tank.hcoef = 150;
  tank.Pa = 1e5;
  tank.Ta = 300;
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
          inlet.m_flow = cv*max(0, charF)*sqrt(max(0, med.d))*sqrt(max(inlet.p, outlet.p))*sign(inlet.p - outlet.p)*regRoot(1 - max(prat, 0.5), dpTol);
        else
          inlet.m_flow = cv*charF*sqrt(med.d*inlet.p)*regRoot(1 - outlet.p/inlet.p, dpTol);
        end if;
      end Valve;

      type Vchar = enumeration(Linear, FastActing, EquiPercent);
    end Valves;

    package Tanks
      model SimpleTank
        replaceable package Medium = PartialMixtureMedium;
        FluidPort inlet(redeclare package Medium = Medium);
        FluidPort outlet(redeclare package Medium = Medium);
        Heat Q_in;
        Temperature Tf;
        Pressure Pa;
        Energy U;
        SpecificHeatCapacity Cp;
        Medium.ThermodynamicState state;
        Mass mf;
        Volume fvol;
        Volume tvol;
        Density rho;
      equation
        state = Medium.setState_pTX(Pa, Tf, outlet.Xi_outflow);
        outlet.h_outflow = Medium.specificEnthalpy(state);
        rho = Medium.density(state);
        der(mf) = smooth(1, if (fvol < tvol) then (inlet.m_flow + outlet.m_flow) else 0);
        fvol = mf/rho;
        inlet.m_flow*actualStream(inlet.h_outflow) + outlet.m_flow*actualStream(outlet.h_outflow) + Q_in = der(U);
        U = fvol*Medium.density(state)*Medium.specificInternalEnergy(state);
        Cp = Medium.specificHeatCapacityCp(state);
        outlet.Xi_outflow = inStream(inlet.Xi_outflow);
        inlet.Xi_outflow = inStream(outlet.Xi_outflow);
        inlet.h_outflow = outlet.h_outflow;
      end SimpleTank;

      model OpenTank
        extends ThermoS.Uops.Interfaces.CylindricalTank;
        extends SimpleTank;
        parameter Length in_pos = h;
        CoefficientOfHeatTransfer hcoef;
        Area A_h;
        Temperature Ta;
      equation
        outlet.p = Pa + rho*9.8*level;
        assert(level <= h or inlet.m_flow > 0, "Tank overflowing  at t=" + String(time), AssertionLevel.warning);
        assert(level > 0, "Tank Empty  at t=" + String(time));
        inlet.p = smooth(1, if level > in_pos then Pa + rho*9.8*(level - in_pos) else Pa);
        assert(level < in_pos, "Tank Level is above the Inlet at t=" + String(time), AssertionLevel.warning);
        assert(level > in_pos, "Tank Level is below the Inlet at t=" + String(time), AssertionLevel.warning);
        assert((level < in_pos) or (level > in_pos), "Tank Level is at the Inlet at t=" + String(time), AssertionLevel.warning);
        A_h = wettedArea;
        Q_in = hcoef*A_h*(Ta - Tf);
      end OpenTank;
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

    package Interfaces
      model CylindricalTank
        import Modelica.Constants.*;
        parameter Length h = 1.0;
        parameter Length d = 1.0;
        Volume fvol;
        Volume tvol;
        Length level;
        Percent pFull;
        Area wettedArea;
        Area interfaceArea;
      equation
        fvol = pi*d*d*level/4;
        tvol = pi*d*d*h/4;
        pFull = 100*fvol/tvol;
        wettedArea = (pi*d*d/4) + (pi*d*level);
        interfaceArea = pi*d*d/4;
      end CylindricalTank;
    end Interfaces;

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

    import Modelica.Media.Interfaces.PartialMixtureMedium;
    import Modelica.Fluid.Interfaces.FluidPort;
    import Modelica.Fluid.Utilities.*;
    import Modelica.Units.SI.*;
    import ThermoS.Types.*;
  end Uops;

  package Types
    type Percent = Real(unit = "%", min = 0, max = 100);
    type Fraction = Real(min = 0, max = 1.0);
  end Types;

  package Media
    package MyGas
      import Modelica.Media.IdealGases.Common.MixtureGasNasa;
      import Modelica.Media.IdealGases.Common.SingleGasesData;
      import Modelica.Media.IdealGases.Common.FluidData;
      extends MixtureGasNasa(data = {SingleGasesData.N2, // note data is of type  DataRecord[:]
      SingleGasesData.O2, SingleGasesData.CO2}, fluidConstants = {FluidData.N2, FluidData.O2, FluidData.CO2}, substanceNames = {"Nitrogen", "Oxygen", "CarbonDioxide"}, reducedX = true, reference_X = {0.7, 0.2, 0.1}, Density(start = 1, nominal = 1), AbsolutePressure(start = 1e5, min = 1e3, max = 50e5, nominal = 1e5), Temperature(start = 300, min = 200, max = 2000, nominal = 300), MassFraction(start = 0.333333), MoleFraction(start = 0.333333));
    end MyGas;
  end Media;
end ThermoS;

model plant_total
  extends plant;
end plant_total;
