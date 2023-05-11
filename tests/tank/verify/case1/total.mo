package ModelicaServices
  constant String target = "OpenModelica";

  package Machine
    final constant Real eps = 1e-15;
    final constant Real small = 1e-60;
    final constant Real inf = 1e60;
    final constant Integer Integer_inf = OpenModelica.Internal.Architecture.integerMax();
  end Machine;

  package Types end Types;
  annotation(version = "4.0.0", versionDate = "2020-06-04", dateModified = "2020-06-04 11:00:00Z");
end ModelicaServices;

package Modelica
  package Fluid
    import Modelica.Units.SI;
    import Cv = Modelica.Units.Conversions;

    package Interfaces
      connector FluidPort
        replaceable package Medium = Modelica.Media.Interfaces.PartialMedium annotation(choicesAllMatching = true);
        flow Medium.MassFlowRate m_flow;
        Medium.AbsolutePressure p;
        stream Medium.SpecificEnthalpy h_outflow;
        stream Medium.MassFraction[Medium.nXi] Xi_outflow;
        stream Medium.ExtraProperty[Medium.nC] C_outflow;
      end FluidPort;
    end Interfaces;

    package Types end Types;

    package Utilities end Utilities;
  end Fluid;

  package Media
    import Modelica.Units.SI;
    import Cv = Modelica.Units.Conversions;

    package Interfaces
      partial package PartialMedium
        extends Modelica.Media.Interfaces.Types;
        constant Modelica.Media.Interfaces.Choices.IndependentVariables ThermoStates;
        constant String mediumName = "unusablePartialMedium";
        constant String[:] substanceNames = {mediumName};
        constant String[:] extraPropertiesNames = fill("", 0);
        constant Boolean singleState;
        constant Boolean reducedX = true;
        constant Boolean fixedX = false;
        constant AbsolutePressure reference_p = 101325;
        constant Temperature reference_T = 298.15;
        constant MassFraction[nX] reference_X = fill(1/nX, nX);
        constant AbsolutePressure p_default = 101325;
        constant Temperature T_default = Modelica.Units.Conversions.from_degC(20);
        constant SpecificEnthalpy h_default = specificEnthalpy_pTX(p_default, T_default, X_default);
        constant MassFraction[nX] X_default = reference_X;
        constant ExtraProperty[nC] C_default = fill(0, nC);
        final constant Integer nS = size(substanceNames, 1);
        constant Integer nX = nS;
        constant Integer nXi = if fixedX then 0 else if reducedX then nS - 1 else nS;
        final constant Integer nC = size(extraPropertiesNames, 1);
        constant Real[nC] C_nominal(min = fill(Modelica.Constants.eps, nC)) = 1.0e-6*ones(nC);
        replaceable record FluidConstants = Modelica.Media.Interfaces.Types.Basic.FluidConstants;

        replaceable record ThermodynamicState end ThermodynamicState;

        replaceable partial model BaseProperties
          InputAbsolutePressure p;
          InputMassFraction[nXi] Xi(start = reference_X[1:nXi]);
          InputSpecificEnthalpy h;
          Density d;
          Temperature T;
          MassFraction[nX] X(start = reference_X);
          SpecificInternalEnergy u;
          SpecificHeatCapacity R_s;
          MolarMass MM;
          ThermodynamicState state;
          parameter Boolean preferredMediumStates = false annotation(Evaluate = true);
          parameter Boolean standardOrderComponents = true;
          Modelica.Units.NonSI.Temperature_degC T_degC = Modelica.Units.Conversions.to_degC(T);
          Modelica.Units.NonSI.Pressure_bar p_bar = Modelica.Units.Conversions.to_bar(p);
          connector InputAbsolutePressure = input SI.AbsolutePressure;
          connector InputSpecificEnthalpy = input SI.SpecificEnthalpy;
          connector InputMassFraction = input SI.MassFraction;
        equation
          if standardOrderComponents then
            Xi = X[1:nXi];
            if fixedX then
              X = reference_X;
            end if;
            if reducedX and not fixedX then
              X[nX] = 1 - sum(Xi);
            end if;
            for i in 1:nX loop
              assert(X[i] >= -1.e-5 and X[i] <= 1 + 1.e-5, "Mass fraction X[" + String(i) + "] = " + String(X[i]) + "of substance " + substanceNames[i] + "\nof medium " + mediumName + " is not in the range 0..1");
            end for;
          end if;
          assert(p >= 0.0, "Pressure (= " + String(p) + " Pa) of medium \"" + mediumName + "\" is negative\n(Temperature = " + String(T) + " K)");
        end BaseProperties;

        replaceable partial function setState_pTX
          input AbsolutePressure p;
          input Temperature T;
          input MassFraction[:] X = reference_X;
          output ThermodynamicState state;
        end setState_pTX;

        replaceable partial function setState_phX
          input AbsolutePressure p;
          input SpecificEnthalpy h;
          input MassFraction[:] X = reference_X;
          output ThermodynamicState state;
        end setState_phX;

        replaceable partial function setState_dTX
          input Density d;
          input Temperature T;
          input MassFraction[:] X = reference_X;
          output ThermodynamicState state;
        end setState_dTX;

        replaceable partial function pressure
          input ThermodynamicState state;
          output AbsolutePressure p;
        end pressure;

        replaceable partial function temperature
          input ThermodynamicState state;
          output Temperature T;
        end temperature;

        replaceable partial function density
          input ThermodynamicState state;
          output Density d;
        end density;

        replaceable partial function specificEnthalpy
          input ThermodynamicState state;
          output SpecificEnthalpy h;
        end specificEnthalpy;

        replaceable partial function specificInternalEnergy
          input ThermodynamicState state;
          output SpecificEnergy u;
        end specificInternalEnergy;

        replaceable partial function specificHeatCapacityCp
          input ThermodynamicState state;
          output SpecificHeatCapacity cp;
        end specificHeatCapacityCp;

        replaceable partial function molarMass
          input ThermodynamicState state;
          output MolarMass MM;
        end molarMass;

        replaceable function specificEnthalpy_pTX
          input AbsolutePressure p;
          input Temperature T;
          input MassFraction[:] X = reference_X;
          output SpecificEnthalpy h;
        algorithm
          h := specificEnthalpy(setState_pTX(p, T, X));
          annotation(inverse(T = temperature_phX(p, h, X)));
        end specificEnthalpy_pTX;

        replaceable function temperature_phX
          input AbsolutePressure p;
          input SpecificEnthalpy h;
          input MassFraction[:] X = reference_X;
          output Temperature T;
        algorithm
          T := temperature(setState_phX(p, h, X));
        end temperature_phX;

        replaceable function density_phX
          input AbsolutePressure p;
          input SpecificEnthalpy h;
          input MassFraction[:] X = reference_X;
          output Density d;
        algorithm
          d := density(setState_phX(p, h, X));
        end density_phX;

        type MassFlowRate = SI.MassFlowRate(quantity = "MassFlowRate." + mediumName, min = -1.0e5, max = 1.e5);
      end PartialMedium;

      partial package PartialPureSubstance
        extends PartialMedium(final reducedX = true, final fixedX = true);

        replaceable function density_ph
          input AbsolutePressure p;
          input SpecificEnthalpy h;
          output Density d;
        algorithm
          d := density_phX(p, h, fill(0, 0));
        end density_ph;

        replaceable function temperature_ph
          input AbsolutePressure p;
          input SpecificEnthalpy h;
          output Temperature T;
        algorithm
          T := temperature_phX(p, h, fill(0, 0));
        end temperature_ph;

        replaceable function pressure_dT
          input Density d;
          input Temperature T;
          output AbsolutePressure p;
        algorithm
          p := pressure(setState_dTX(d, T, fill(0, 0)));
        end pressure_dT;

        replaceable function specificEnthalpy_dT
          input Density d;
          input Temperature T;
          output SpecificEnthalpy h;
        algorithm
          h := specificEnthalpy(setState_dTX(d, T, fill(0, 0)));
        end specificEnthalpy_dT;

        replaceable function specificEnthalpy_pT
          input AbsolutePressure p;
          input Temperature T;
          output SpecificEnthalpy h;
        algorithm
          h := specificEnthalpy_pTX(p, T, fill(0, 0));
        end specificEnthalpy_pT;

        replaceable function density_pT
          input AbsolutePressure p;
          input Temperature T;
          output Density d;
        algorithm
          d := density(setState_pTX(p, T, fill(0, 0)));
        end density_pT;

        redeclare replaceable partial model extends BaseProperties(final standardOrderComponents = true) end BaseProperties;
      end PartialPureSubstance;

      partial package PartialMixtureMedium
        extends PartialMedium(redeclare replaceable record FluidConstants = Modelica.Media.Interfaces.Types.IdealGas.FluidConstants);

        redeclare replaceable record extends ThermodynamicState
          AbsolutePressure p;
          Temperature T;
          MassFraction[nX] X(start = reference_X);
        end ThermodynamicState;

        constant FluidConstants[nS] fluidConstants;
      end PartialMixtureMedium;

      partial package PartialTwoPhaseMedium
        extends PartialPureSubstance(redeclare replaceable record FluidConstants = Modelica.Media.Interfaces.Types.TwoPhase.FluidConstants);
        constant Boolean smoothModel = false;
        constant Boolean onePhase = false;
        constant FluidConstants[nS] fluidConstants;

        redeclare replaceable record extends ThermodynamicState
          FixedPhase phase(min = 0, max = 2);
        end ThermodynamicState;

        redeclare replaceable partial model extends BaseProperties
          SaturationProperties sat;
        end BaseProperties;

        redeclare replaceable partial function extends setState_dTX
          input FixedPhase phase = 0;
        end setState_dTX;

        redeclare replaceable partial function extends setState_phX
          input FixedPhase phase = 0;
        end setState_phX;

        redeclare replaceable partial function extends setState_pTX
          input FixedPhase phase = 0;
        end setState_pTX;

        replaceable partial function bubbleEnthalpy
          input SaturationProperties sat;
          output SI.SpecificEnthalpy hl;
        end bubbleEnthalpy;

        replaceable partial function dewEnthalpy
          input SaturationProperties sat;
          output SI.SpecificEnthalpy hv;
        end dewEnthalpy;

        replaceable partial function bubbleDensity
          input SaturationProperties sat;
          output Density dl;
        end bubbleDensity;

        replaceable partial function dewDensity
          input SaturationProperties sat;
          output Density dv;
        end dewDensity;

        replaceable partial function saturationPressure
          input Temperature T;
          output AbsolutePressure p;
        end saturationPressure;

        replaceable partial function saturationTemperature
          input AbsolutePressure p;
          output Temperature T;
        end saturationTemperature;

        redeclare replaceable function extends molarMass
        algorithm
          MM := fluidConstants[1].molarMass;
        end molarMass;

        redeclare replaceable function specificEnthalpy_pTX
          input AbsolutePressure p;
          input Temperature T;
          input MassFraction[:] X;
          input FixedPhase phase = 0;
          output SpecificEnthalpy h;
        algorithm
          h := specificEnthalpy(setState_pTX(p, T, X, phase));
        end specificEnthalpy_pTX;

        redeclare replaceable function temperature_phX
          input AbsolutePressure p;
          input SpecificEnthalpy h;
          input MassFraction[:] X;
          input FixedPhase phase = 0;
          output Temperature T;
        algorithm
          T := temperature(setState_phX(p, h, X, phase));
        end temperature_phX;

        redeclare replaceable function density_phX
          input AbsolutePressure p;
          input SpecificEnthalpy h;
          input MassFraction[:] X;
          input FixedPhase phase = 0;
          output Density d;
        algorithm
          d := density(setState_phX(p, h, X, phase));
        end density_phX;

        redeclare replaceable function density_ph
          input AbsolutePressure p;
          input SpecificEnthalpy h;
          input FixedPhase phase = 0;
          output Density d;
        algorithm
          d := density_phX(p, h, fill(0, 0), phase);
        end density_ph;

        redeclare replaceable function temperature_ph
          input AbsolutePressure p;
          input SpecificEnthalpy h;
          input FixedPhase phase = 0;
          output Temperature T;
        algorithm
          T := temperature_phX(p, h, fill(0, 0), phase);
        end temperature_ph;

        redeclare replaceable function pressure_dT
          input Density d;
          input Temperature T;
          input FixedPhase phase = 0;
          output AbsolutePressure p;
        algorithm
          p := pressure(setState_dTX(d, T, fill(0, 0), phase));
        end pressure_dT;

        redeclare replaceable function specificEnthalpy_dT
          input Density d;
          input Temperature T;
          input FixedPhase phase = 0;
          output SpecificEnthalpy h;
        algorithm
          h := specificEnthalpy(setState_dTX(d, T, fill(0, 0), phase));
        end specificEnthalpy_dT;

        redeclare replaceable function specificEnthalpy_pT
          input AbsolutePressure p;
          input Temperature T;
          input FixedPhase phase = 0;
          output SpecificEnthalpy h;
        algorithm
          h := specificEnthalpy_pTX(p, T, fill(0, 0), phase);
        end specificEnthalpy_pT;

        redeclare replaceable function density_pT
          input AbsolutePressure p;
          input Temperature T;
          input FixedPhase phase = 0;
          output Density d;
        algorithm
          d := density(setState_pTX(p, T, fill(0, 0), phase));
        end density_pT;
      end PartialTwoPhaseMedium;

      partial package PartialSimpleMedium
        extends Interfaces.PartialPureSubstance(final ThermoStates = Modelica.Media.Interfaces.Choices.IndependentVariables.pT, final singleState = true);
        constant SpecificHeatCapacity cp_const;
        constant SpecificHeatCapacity cv_const;
        constant Density d_const;
        constant DynamicViscosity eta_const;
        constant ThermalConductivity lambda_const;
        constant VelocityOfSound a_const;
        constant Temperature T_min;
        constant Temperature T_max;
        constant Temperature T0 = reference_T;
        constant MolarMass MM_const;
        constant FluidConstants[nS] fluidConstants;

        redeclare record extends ThermodynamicState
          AbsolutePressure p;
          Temperature T;
        end ThermodynamicState;

        redeclare replaceable model extends BaseProperties(T(stateSelect = if preferredMediumStates then StateSelect.prefer else StateSelect.default), p(stateSelect = if preferredMediumStates then StateSelect.prefer else StateSelect.default))
        equation
          assert(T >= T_min and T <= T_max, "
        Temperature T (= " + String(T) + " K) is not
        in the allowed range (" + String(T_min) + " K <= T <= " + String(T_max) + " K)
        required from medium model \"" + mediumName + "\".
          ");
          h = specificEnthalpy_pTX(p, T, X);
          u = cv_const*(T - T0);
          d = d_const;
          R_s = 0;
          MM = MM_const;
          state.T = T;
          state.p = p;
        end BaseProperties;

        redeclare function setState_pTX
          input AbsolutePressure p;
          input Temperature T;
          input MassFraction[:] X = reference_X;
          output ThermodynamicState state;
        algorithm
          state := ThermodynamicState(p = p, T = T);
        end setState_pTX;

        redeclare function setState_phX
          input AbsolutePressure p;
          input SpecificEnthalpy h;
          input MassFraction[:] X = reference_X;
          output ThermodynamicState state;
        algorithm
          state := ThermodynamicState(p = p, T = T0 + h/cp_const);
        end setState_phX;

        redeclare function setState_dTX
          input Density d;
          input Temperature T;
          input MassFraction[:] X = reference_X;
          output ThermodynamicState state;
        algorithm
          assert(false, "Pressure can not be computed from temperature and density for an incompressible fluid!");
        end setState_dTX;

        redeclare function extends pressure
        algorithm
          p := state.p;
        end pressure;

        redeclare function extends temperature
        algorithm
          T := state.T;
        end temperature;

        redeclare function extends density
        algorithm
          d := d_const;
        end density;

        redeclare function extends specificEnthalpy
        algorithm
          h := cp_const*(state.T - T0);
        end specificEnthalpy;

        redeclare function extends specificHeatCapacityCp
        algorithm
          cp := cp_const;
        end specificHeatCapacityCp;

        redeclare function specificEnthalpy_pTX
          input AbsolutePressure p;
          input Temperature T;
          input MassFraction[nX] X;
          output SpecificEnthalpy h;
        algorithm
          h := cp_const*(T - T0);
        end specificEnthalpy_pTX;

        redeclare function temperature_phX
          input AbsolutePressure p;
          input SpecificEnthalpy h;
          input MassFraction[nX] X;
          output Temperature T;
        algorithm
          T := T0 + h/cp_const;
        end temperature_phX;

        redeclare function density_phX
          input AbsolutePressure p;
          input SpecificEnthalpy h;
          input MassFraction[nX] X;
          output Density d;
        algorithm
          d := density(setState_phX(p, h, X));
        end density_phX;

        redeclare function extends specificInternalEnergy
        algorithm
          u := cv_const*(state.T - T0);
        end specificInternalEnergy;

        redeclare function extends molarMass
        algorithm
          MM := MM_const;
        end molarMass;
      end PartialSimpleMedium;

      package Choices
        type IndependentVariables = enumeration(T, pT, ph, phX, pTX, dTX);
        type pd = enumeration(default, p_known, d_known) annotation(Evaluate = true);
      end Choices;

      package Types
        type AbsolutePressure = SI.AbsolutePressure(min = 0, max = 1.e8, nominal = 1.e5, start = 1.e5);
        type Density = SI.Density(min = 0, max = 1.e5, nominal = 1, start = 1);
        type DynamicViscosity = SI.DynamicViscosity(min = 0, max = 1.e8, nominal = 1.e-3, start = 1.e-3);
        type MassFraction = Real(quantity = "MassFraction", final unit = "kg/kg", min = 0, max = 1, nominal = 0.1);
        type MoleFraction = Real(quantity = "MoleFraction", final unit = "mol/mol", min = 0, max = 1, nominal = 0.1);
        type MolarMass = SI.MolarMass(min = 0.001, max = 0.25, nominal = 0.032);
        type MolarVolume = SI.MolarVolume(min = 1e-6, max = 1.0e6, nominal = 1.0);
        type SpecificEnergy = SI.SpecificEnergy(min = -1.0e8, max = 1.e8, nominal = 1.e6);
        type SpecificInternalEnergy = SpecificEnergy;
        type SpecificEnthalpy = SI.SpecificEnthalpy(min = -1.0e10, max = 1.e10, nominal = 1.e6);
        type SpecificEntropy = SI.SpecificEntropy(min = -1.e7, max = 1.e7, nominal = 1.e3);
        type SpecificHeatCapacity = SI.SpecificHeatCapacity(min = 0, max = 1.e7, nominal = 1.e3, start = 1.e3);
        type Temperature = SI.Temperature(min = 1, max = 1.e4, nominal = 300, start = 288.15);
        type ThermalConductivity = SI.ThermalConductivity(min = 0, max = 500, nominal = 1, start = 1);
        type VelocityOfSound = SI.Velocity(min = 0, max = 1.e5, nominal = 1000, start = 1000);
        type ExtraProperty = Real(min = 0.0, start = 1.0);
        type DipoleMoment = Real(min = 0.0, max = 2.0, unit = "debye", quantity = "ElectricDipoleMoment");

        replaceable record SaturationProperties
          AbsolutePressure psat;
          Temperature Tsat;
        end SaturationProperties;

        type FixedPhase = Integer(min = 0, max = 2);

        package Basic
          record FluidConstants
            String iupacName;
            String casRegistryNumber;
            String chemicalFormula;
            String structureFormula;
            MolarMass molarMass;
          end FluidConstants;
        end Basic;

        package IdealGas
          record FluidConstants
            extends Modelica.Media.Interfaces.Types.Basic.FluidConstants;
            Temperature criticalTemperature;
            AbsolutePressure criticalPressure;
            MolarVolume criticalMolarVolume;
            Real acentricFactor;
            Temperature meltingPoint;
            Temperature normalBoilingPoint;
            DipoleMoment dipoleMoment;
            Boolean hasIdealGasHeatCapacity = false;
            Boolean hasCriticalData = false;
            Boolean hasDipoleMoment = false;
            Boolean hasFundamentalEquation = false;
            Boolean hasLiquidHeatCapacity = false;
            Boolean hasSolidHeatCapacity = false;
            Boolean hasAccurateViscosityData = false;
            Boolean hasAccurateConductivityData = false;
            Boolean hasVapourPressureCurve = false;
            Boolean hasAcentricFactor = false;
            SpecificEnthalpy HCRIT0 = 0.0;
            SpecificEntropy SCRIT0 = 0.0;
            SpecificEnthalpy deltah = 0.0;
            SpecificEntropy deltas = 0.0;
          end FluidConstants;
        end IdealGas;

        package TwoPhase
          record FluidConstants
            extends Modelica.Media.Interfaces.Types.Basic.FluidConstants;
            Temperature criticalTemperature;
            AbsolutePressure criticalPressure;
            MolarVolume criticalMolarVolume;
            Real acentricFactor;
            Temperature triplePointTemperature;
            AbsolutePressure triplePointPressure;
            Temperature meltingPoint;
            Temperature normalBoilingPoint;
            DipoleMoment dipoleMoment;
            Boolean hasIdealGasHeatCapacity = false;
            Boolean hasCriticalData = false;
            Boolean hasDipoleMoment = false;
            Boolean hasFundamentalEquation = false;
            Boolean hasLiquidHeatCapacity = false;
            Boolean hasSolidHeatCapacity = false;
            Boolean hasAccurateViscosityData = false;
            Boolean hasAccurateConductivityData = false;
            Boolean hasVapourPressureCurve = false;
            Boolean hasAcentricFactor = false;
            SpecificEnthalpy HCRIT0 = 0.0;
            SpecificEntropy SCRIT0 = 0.0;
            SpecificEnthalpy deltah = 0.0;
            SpecificEntropy deltas = 0.0;
          end FluidConstants;
        end TwoPhase;
      end Types;
    end Interfaces;

    package Common
      type DerPressureByDensity = Real(final quantity = "DerPressureByDensity", final unit = "Pa.m3/kg");
      type DerPressureByTemperature = Real(final quantity = "DerPressureByTemperature", final unit = "Pa/K");
      constant Real MINPOS = 1.0e-9;
      constant SI.Area AMIN = MINPOS;
      constant SI.Area AMAX = 1.0e5;
      constant SI.Area ANOM = 1.0;
      constant SI.AmountOfSubstance MOLMIN = -1.0*MINPOS;
      constant SI.AmountOfSubstance MOLMAX = 1.0e8;
      constant SI.AmountOfSubstance MOLNOM = 1.0;
      constant SI.Density DMIN = 1e-6;
      constant SI.Density DMAX = 30.0e3;
      constant SI.Density DNOM = 1.0;
      constant SI.ThermalConductivity LAMMIN = MINPOS;
      constant SI.ThermalConductivity LAMNOM = 1.0;
      constant SI.ThermalConductivity LAMMAX = 1000.0;
      constant SI.DynamicViscosity ETAMIN = MINPOS;
      constant SI.DynamicViscosity ETAMAX = 1.0e8;
      constant SI.DynamicViscosity ETANOM = 100.0;
      constant SI.Energy EMIN = -1.0e10;
      constant SI.Energy EMAX = 1.0e10;
      constant SI.Energy ENOM = 1.0e3;
      constant SI.Entropy SMIN = -1.0e6;
      constant SI.Entropy SMAX = 1.0e6;
      constant SI.Entropy SNOM = 1.0e3;
      constant SI.MassFlowRate MDOTMIN = -1.0e5;
      constant SI.MassFlowRate MDOTMAX = 1.0e5;
      constant SI.MassFlowRate MDOTNOM = 1.0;
      constant SI.MassFraction MASSXMIN = -1.0*MINPOS;
      constant SI.MassFraction MASSXMAX = 1.0;
      constant SI.MassFraction MASSXNOM = 0.1;
      constant SI.Mass MMIN = -1.0*MINPOS;
      constant SI.Mass MMAX = 1.0e8;
      constant SI.Mass MNOM = 1.0;
      constant SI.MolarMass MMMIN = 0.001;
      constant SI.MolarMass MMMAX = 250.0;
      constant SI.MolarMass MMNOM = 0.2;
      constant SI.MoleFraction MOLEYMIN = -1.0*MINPOS;
      constant SI.MoleFraction MOLEYMAX = 1.0;
      constant SI.MoleFraction MOLEYNOM = 0.1;
      constant SI.MomentumFlux GMIN = -1.0e8;
      constant SI.MomentumFlux GMAX = 1.0e8;
      constant SI.MomentumFlux GNOM = 1.0;
      constant SI.Power POWMIN = -1.0e8;
      constant SI.Power POWMAX = 1.0e8;
      constant SI.Power POWNOM = 1.0e3;
      constant SI.Pressure PMIN = 1.0e4;
      constant SI.Pressure PMAX = 1.0e8;
      constant SI.Pressure PNOM = 1.0e5;
      constant SI.Pressure COMPPMIN = -1.0*MINPOS;
      constant SI.Pressure COMPPMAX = 1.0e8;
      constant SI.Pressure COMPPNOM = 1.0e5;
      constant SI.RatioOfSpecificHeatCapacities KAPPAMIN = 1.0;
      constant SI.RatioOfSpecificHeatCapacities KAPPAMAX = 1.7;
      constant SI.RatioOfSpecificHeatCapacities KAPPANOM = 1.2;
      constant SI.SpecificEnergy SEMIN = -1.0e8;
      constant SI.SpecificEnergy SEMAX = 1.0e8;
      constant SI.SpecificEnergy SENOM = 1.0e6;
      constant SI.SpecificEnthalpy SHMIN = -1.0e8;
      constant SI.SpecificEnthalpy SHMAX = 1.0e8;
      constant SI.SpecificEnthalpy SHNOM = 1.0e6;
      constant SI.SpecificEntropy SSMIN = -1.0e6;
      constant SI.SpecificEntropy SSMAX = 1.0e6;
      constant SI.SpecificEntropy SSNOM = 1.0e3;
      constant SI.SpecificHeatCapacity CPMIN = MINPOS;
      constant SI.SpecificHeatCapacity CPMAX = 1.0e6;
      constant SI.SpecificHeatCapacity CPNOM = 1.0e3;
      constant SI.Temperature TMIN = 1.0;
      constant SI.Temperature TMAX = 6000.0;
      constant SI.Temperature TNOM = 320.0;
      constant SI.ThermalConductivity LMIN = MINPOS;
      constant SI.ThermalConductivity LMAX = 500.0;
      constant SI.ThermalConductivity LNOM = 1.0;
      constant SI.Velocity VELMIN = -1.0e5;
      constant SI.Velocity VELMAX = 1.0e5;
      constant SI.Velocity VELNOM = 1.0;
      constant SI.Volume VMIN = 0.0;
      constant SI.Volume VMAX = 1.0e5;
      constant SI.Volume VNOM = 1.0e-3;

      record SaturationProperties
        SI.Temperature T;
        SI.Density d;
        SI.Pressure p;
        SI.SpecificEnergy u;
        SI.SpecificEnthalpy h;
        SI.SpecificEntropy s;
        SI.SpecificHeatCapacity cp;
        SI.SpecificHeatCapacity cv;
        SI.SpecificHeatCapacity R_s;
        SI.RatioOfSpecificHeatCapacities kappa;
        PhaseBoundaryProperties liq;
        PhaseBoundaryProperties vap;
        Real dpT(unit = "Pa/K");
        SI.MassFraction x;
      end SaturationProperties;

      record IF97BaseTwoPhase
        Integer phase(start = 0);
        Integer region(min = 1, max = 5);
        SI.Pressure p;
        SI.Temperature T;
        SI.SpecificEnthalpy h;
        SI.SpecificHeatCapacity R_s;
        SI.SpecificHeatCapacity cp;
        SI.SpecificHeatCapacity cv;
        SI.Density rho;
        SI.SpecificEntropy s;
        DerPressureByTemperature pt;
        DerPressureByDensity pd;
        Real vt;
        Real vp;
        Real x;
        Real dpT;
      end IF97BaseTwoPhase;

      record IF97PhaseBoundaryProperties
        Boolean region3boundary;
        SI.SpecificHeatCapacity R_s;
        SI.Temperature T;
        SI.Density d;
        SI.SpecificEnthalpy h;
        SI.SpecificEntropy s;
        SI.SpecificHeatCapacity cp;
        SI.SpecificHeatCapacity cv;
        DerPressureByTemperature dpT;
        DerPressureByTemperature pt;
        DerPressureByDensity pd;
        Real vt(unit = "m3/(kg.K)");
        Real vp(unit = "m3/(kg.Pa)");
      end IF97PhaseBoundaryProperties;

      record GibbsDerivs
        SI.Pressure p;
        SI.Temperature T;
        SI.SpecificHeatCapacity R_s;
        Real pi(unit = "1");
        Real tau(unit = "1");
        Real g(unit = "1");
        Real gpi(unit = "1");
        Real gpipi(unit = "1");
        Real gtau(unit = "1");
        Real gtautau(unit = "1");
        Real gtaupi(unit = "1");
      end GibbsDerivs;

      record HelmholtzDerivs
        SI.Density d;
        SI.Temperature T;
        SI.SpecificHeatCapacity R_s;
        Real delta(unit = "1");
        Real tau(unit = "1");
        Real f(unit = "1");
        Real fdelta(unit = "1");
        Real fdeltadelta(unit = "1");
        Real ftau(unit = "1");
        Real ftautau(unit = "1");
        Real fdeltatau(unit = "1");
      end HelmholtzDerivs;

      record PhaseBoundaryProperties
        SI.Density d;
        SI.SpecificEnthalpy h;
        SI.SpecificEnergy u;
        SI.SpecificEntropy s;
        SI.SpecificHeatCapacity cp;
        SI.SpecificHeatCapacity cv;
        DerPressureByTemperature pt;
        DerPressureByDensity pd;
      end PhaseBoundaryProperties;

      record NewtonDerivatives_ph
        SI.Pressure p;
        SI.SpecificEnthalpy h;
        DerPressureByDensity pd;
        DerPressureByTemperature pt;
        Real hd;
        Real ht;
      end NewtonDerivatives_ph;

      record NewtonDerivatives_pT
        SI.Pressure p;
        DerPressureByDensity pd;
      end NewtonDerivatives_pT;

      function gibbsToBoundaryProps
        input GibbsDerivs g;
        output PhaseBoundaryProperties sat;
      protected
        Real vt;
        Real vp;
      algorithm
        sat.d := g.p/(g.R_s*g.T*g.pi*g.gpi);
        sat.h := g.R_s*g.T*g.tau*g.gtau;
        sat.u := g.T*g.R_s*(g.tau*g.gtau - g.pi*g.gpi);
        sat.s := g.R_s*(g.tau*g.gtau - g.g);
        sat.cp := -g.R_s*g.tau*g.tau*g.gtautau;
        sat.cv := g.R_s*(-g.tau*g.tau*g.gtautau + (g.gpi - g.tau*g.gtaupi)*(g.gpi - g.tau*g.gtaupi)/(g.gpipi));
        vt := g.R_s/g.p*(g.pi*g.gpi - g.tau*g.pi*g.gtaupi);
        vp := g.R_s*g.T/(g.p*g.p)*g.pi*g.pi*g.gpipi;
        sat.pt := -g.p/g.T*(g.gpi - g.tau*g.gtaupi)/(g.gpipi*g.pi);
        sat.pd := -g.R_s*g.T*g.gpi*g.gpi/(g.gpipi);
      end gibbsToBoundaryProps;

      function helmholtzToBoundaryProps
        input HelmholtzDerivs f;
        output PhaseBoundaryProperties sat;
      protected
        SI.Pressure p;
      algorithm
        p := f.R_s*f.d*f.T*f.delta*f.fdelta;
        sat.d := f.d;
        sat.h := f.R_s*f.T*(f.tau*f.ftau + f.delta*f.fdelta);
        sat.s := f.R_s*(f.tau*f.ftau - f.f);
        sat.u := f.R_s*f.T*f.tau*f.ftau;
        sat.cp := f.R_s*(-f.tau*f.tau*f.ftautau + (f.delta*f.fdelta - f.delta*f.tau*f.fdeltatau)^2/(2*f.delta*f.fdelta + f.delta*f.delta*f.fdeltadelta));
        sat.cv := f.R_s*(-f.tau*f.tau*f.ftautau);
        sat.pt := f.R_s*f.d*f.delta*(f.fdelta - f.tau*f.fdeltatau);
        sat.pd := f.R_s*f.T*f.delta*(2.0*f.fdelta + f.delta*f.fdeltadelta);
      end helmholtzToBoundaryProps;

      function cv2Phase
        input PhaseBoundaryProperties liq;
        input PhaseBoundaryProperties vap;
        input SI.MassFraction x;
        input SI.Temperature T;
        input SI.Pressure p;
        output SI.SpecificHeatCapacity cv;
      protected
        Real dpT;
        Real dxv;
        Real dvTl;
        Real dvTv;
        Real duTl;
        Real duTv;
        Real dxt;
      algorithm
        dxv := if (liq.d <> vap.d) then liq.d*vap.d/(liq.d - vap.d) else 0.0;
        dpT := (vap.s - liq.s)*dxv;
        dvTl := (liq.pt - dpT)/liq.pd/liq.d/liq.d;
        dvTv := (vap.pt - dpT)/vap.pd/vap.d/vap.d;
        dxt := -dxv*(dvTl + x*(dvTv - dvTl));
        duTl := liq.cv + (T*liq.pt - p)*dvTl;
        duTv := vap.cv + (T*vap.pt - p)*dvTv;
        cv := duTl + x*(duTv - duTl) + dxt*(vap.u - liq.u);
      end cv2Phase;

      function Helmholtz_ph
        input HelmholtzDerivs f;
        output NewtonDerivatives_ph nderivs;
      protected
        SI.SpecificHeatCapacity cv;
      algorithm
        cv := -f.R_s*(f.tau*f.tau*f.ftautau);
        nderivs.p := f.d*f.R_s*f.T*f.delta*f.fdelta;
        nderivs.h := f.R_s*f.T*(f.tau*f.ftau + f.delta*f.fdelta);
        nderivs.pd := f.R_s*f.T*f.delta*(2.0*f.fdelta + f.delta*f.fdeltadelta);
        nderivs.pt := f.R_s*f.d*f.delta*(f.fdelta - f.tau*f.fdeltatau);
        nderivs.ht := cv + nderivs.pt/f.d;
        nderivs.hd := (nderivs.pd - f.T*nderivs.pt/f.d)/f.d;
      end Helmholtz_ph;

      function Helmholtz_pT
        input HelmholtzDerivs f;
        output NewtonDerivatives_pT nderivs;
      algorithm
        nderivs.p := f.d*f.R_s*f.T*f.delta*f.fdelta;
        nderivs.pd := f.R_s*f.T*f.delta*(2.0*f.fdelta + f.delta*f.fdeltadelta);
      end Helmholtz_pT;
    end Common;

    package Water
      import Modelica.Media.Water.ConstantPropertyLiquidWater.simpleWaterConstants;
      constant Modelica.Media.Interfaces.Types.TwoPhase.FluidConstants[1] waterConstants(each chemicalFormula = "H2O", each structureFormula = "H2O", each casRegistryNumber = "7732-18-5", each iupacName = "oxidane", each molarMass = 0.018015268, each criticalTemperature = 647.096, each criticalPressure = 22064.0e3, each criticalMolarVolume = 1/322.0*0.018015268, each normalBoilingPoint = 373.124, each meltingPoint = 273.15, each triplePointTemperature = 273.16, each triplePointPressure = 611.657, each acentricFactor = 0.344, each dipoleMoment = 1.8, each hasCriticalData = true);

      package ConstantPropertyLiquidWater
        constant Modelica.Media.Interfaces.Types.Basic.FluidConstants[1] simpleWaterConstants(each chemicalFormula = "H2O", each structureFormula = "H2O", each casRegistryNumber = "7732-18-5", each iupacName = "oxidane", each molarMass = 0.018015268);
        extends Interfaces.PartialSimpleMedium(mediumName = "SimpleLiquidWater", cp_const = 4184, cv_const = 4184, d_const = 995.586, eta_const = 1.e-3, lambda_const = 0.598, a_const = 1484, T_min = Cv.from_degC(-1), T_max = Cv.from_degC(130), T0 = 273.15, MM_const = 0.018015268, fluidConstants = simpleWaterConstants);
      end ConstantPropertyLiquidWater;

      package StandardWater = WaterIF97_ph;

      package WaterIF97_ph
        extends WaterIF97_base(ThermoStates = Modelica.Media.Interfaces.Choices.IndependentVariables.ph, final ph_explicit = true, final dT_explicit = false, final pT_explicit = false, smoothModel = false, onePhase = false);
      end WaterIF97_ph;

      partial package WaterIF97_base
        extends Interfaces.PartialTwoPhaseMedium(mediumName = "WaterIF97", substanceNames = {"water"}, singleState = false, SpecificEnthalpy(start = 1.0e5, nominal = 5.0e5), Density(start = 150, nominal = 500), AbsolutePressure(start = 50e5, nominal = 10e5, min = 611.657, max = 100e6), Temperature(start = 500, nominal = 500, min = 273.15, max = 2273.15), smoothModel = false, onePhase = false, fluidConstants = waterConstants);

        redeclare record extends SaturationProperties end SaturationProperties;

        redeclare record extends ThermodynamicState
          SpecificEnthalpy h;
          Density d;
          Temperature T;
          AbsolutePressure p;
        end ThermodynamicState;

        constant Integer Region = 0;
        constant Boolean ph_explicit;
        constant Boolean dT_explicit;
        constant Boolean pT_explicit;

        redeclare replaceable model extends BaseProperties(h(stateSelect = if ph_explicit and preferredMediumStates then StateSelect.prefer else StateSelect.default), d(stateSelect = if dT_explicit and preferredMediumStates then StateSelect.prefer else StateSelect.default), T(stateSelect = if (pT_explicit or dT_explicit) and preferredMediumStates then StateSelect.prefer else StateSelect.default), p(stateSelect = if (pT_explicit or ph_explicit) and preferredMediumStates then StateSelect.prefer else StateSelect.default))
          Integer phase(min = 0, max = 2, start = 1, fixed = false);
        equation
          MM = fluidConstants[1].molarMass;
          if Region > 0 then
            phase = (if Region == 4 then 2 else 1);
          elseif smoothModel then
            if onePhase then
              phase = 1;
              if ph_explicit then
                assert(((h < bubbleEnthalpy(sat) or h > dewEnthalpy(sat)) or p > fluidConstants[1].criticalPressure), "With onePhase=true this model may only be called with one-phase states h < hl or h > hv!" + "(p = " + String(p) + ", h = " + String(h) + ")");
              else
                if dT_explicit then
                  assert(not ((d < bubbleDensity(sat) and d > dewDensity(sat)) and T < fluidConstants[1].criticalTemperature), "With onePhase=true this model may only be called with one-phase states d > dl or d < dv!" + "(d = " + String(d) + ", T = " + String(T) + ")");
                end if;
              end if;
            else
              phase = 0;
            end if;
          else
            if ph_explicit then
              phase = if ((h < bubbleEnthalpy(sat) or h > dewEnthalpy(sat)) or p > fluidConstants[1].criticalPressure) then 1 else 2;
            elseif dT_explicit then
              phase = if not ((d < bubbleDensity(sat) and d > dewDensity(sat)) and T < fluidConstants[1].criticalTemperature) then 1 else 2;
            else
              phase = 1;
            end if;
          end if;
          if dT_explicit then
            p = pressure_dT(d, T, phase, Region);
            h = specificEnthalpy_dT(d, T, phase, Region);
            sat.Tsat = T;
            sat.psat = saturationPressure(T);
          elseif ph_explicit then
            d = density_ph(p, h, phase, Region);
            T = temperature_ph(p, h, phase, Region);
            sat.Tsat = saturationTemperature(p);
            sat.psat = p;
          else
            h = specificEnthalpy_pT(p, T, Region);
            d = density_pT(p, T, Region);
            sat.psat = p;
            sat.Tsat = saturationTemperature(p);
          end if;
          u = h - p/d;
          R_s = Modelica.Constants.R/fluidConstants[1].molarMass;
          h = state.h;
          p = state.p;
          T = state.T;
          d = state.d;
          phase = state.phase;
        end BaseProperties;

        redeclare function density_ph
          input AbsolutePressure p;
          input SpecificEnthalpy h;
          input FixedPhase phase = 0;
          input Integer region = Region;
          output Density d;
        algorithm
          d := IF97_Utilities.rho_ph(p, h, phase, region);
          annotation(Inline = true);
        end density_ph;

        redeclare function temperature_ph
          input AbsolutePressure p;
          input SpecificEnthalpy h;
          input FixedPhase phase = 0;
          input Integer region = Region;
          output Temperature T;
        algorithm
          T := IF97_Utilities.T_ph(p, h, phase, region);
          annotation(Inline = true);
        end temperature_ph;

        redeclare function pressure_dT
          input Density d;
          input Temperature T;
          input FixedPhase phase = 0;
          input Integer region = Region;
          output AbsolutePressure p;
        algorithm
          p := IF97_Utilities.p_dT(d, T, phase, region);
          annotation(Inline = true);
        end pressure_dT;

        redeclare function specificEnthalpy_dT
          input Density d;
          input Temperature T;
          input FixedPhase phase = 0;
          input Integer region = Region;
          output SpecificEnthalpy h;
        algorithm
          h := IF97_Utilities.h_dT(d, T, phase, region);
          annotation(Inline = true);
        end specificEnthalpy_dT;

        redeclare function specificEnthalpy_pT
          input AbsolutePressure p;
          input Temperature T;
          input FixedPhase phase = 0;
          input Integer region = Region;
          output SpecificEnthalpy h;
        algorithm
          h := IF97_Utilities.h_pT(p, T, region);
          annotation(Inline = true);
        end specificEnthalpy_pT;

        redeclare function density_pT
          input AbsolutePressure p;
          input Temperature T;
          input FixedPhase phase = 0;
          input Integer region = Region;
          output Density d;
        algorithm
          d := IF97_Utilities.rho_pT(p, T, region);
          annotation(Inline = true);
        end density_pT;

        redeclare function extends pressure
        algorithm
          p := state.p;
          annotation(Inline = true);
        end pressure;

        redeclare function extends temperature
        algorithm
          T := state.T;
          annotation(Inline = true);
        end temperature;

        redeclare function extends density
        algorithm
          d := state.d;
          annotation(Inline = true);
        end density;

        redeclare function extends specificEnthalpy
        algorithm
          h := state.h;
          annotation(Inline = true);
        end specificEnthalpy;

        redeclare function extends specificInternalEnergy
        algorithm
          u := state.h - state.p/state.d;
          annotation(Inline = true);
        end specificInternalEnergy;

        redeclare function extends specificHeatCapacityCp
        algorithm
          cp := if dT_explicit then IF97_Utilities.cp_dT(state.d, state.T, Region) else if pT_explicit then IF97_Utilities.cp_pT(state.p, state.T, Region) else IF97_Utilities.cp_ph(state.p, state.h, Region);
          annotation(Inline = true);
        end specificHeatCapacityCp;

        redeclare function extends bubbleEnthalpy
        algorithm
          hl := IF97_Utilities.BaseIF97.Regions.hl_p(sat.psat);
          annotation(Inline = true);
        end bubbleEnthalpy;

        redeclare function extends dewEnthalpy
        algorithm
          hv := IF97_Utilities.BaseIF97.Regions.hv_p(sat.psat);
          annotation(Inline = true);
        end dewEnthalpy;

        redeclare function extends bubbleDensity
        algorithm
          dl := if ph_explicit or pT_explicit then IF97_Utilities.BaseIF97.Regions.rhol_p(sat.psat) else IF97_Utilities.BaseIF97.Regions.rhol_T(sat.Tsat);
          annotation(Inline = true);
        end bubbleDensity;

        redeclare function extends dewDensity
        algorithm
          dv := if ph_explicit or pT_explicit then IF97_Utilities.BaseIF97.Regions.rhov_p(sat.psat) else IF97_Utilities.BaseIF97.Regions.rhov_T(sat.Tsat);
          annotation(Inline = true);
        end dewDensity;

        redeclare function extends saturationTemperature
        algorithm
          T := IF97_Utilities.BaseIF97.Basic.tsat(p);
          annotation(Inline = true);
        end saturationTemperature;

        redeclare function extends saturationPressure
        algorithm
          p := IF97_Utilities.BaseIF97.Basic.psat(T);
          annotation(Inline = true);
        end saturationPressure;

        redeclare function extends setState_dTX
          input Integer region = Region;
        algorithm
          state := ThermodynamicState(d = d, T = T, phase = if region == 0 then 0 else if region == 4 then 2 else 1, h = specificEnthalpy_dT(d, T, region = region), p = pressure_dT(d, T, region = region));
          annotation(Inline = true);
        end setState_dTX;

        redeclare function extends setState_phX
          input Integer region = Region;
        algorithm
          state := ThermodynamicState(d = density_ph(p, h, region = region), T = temperature_ph(p, h, region = region), phase = if region == 0 then 0 else if region == 4 then 2 else 1, h = h, p = p);
          annotation(Inline = true);
        end setState_phX;

        redeclare function extends setState_pTX
          input Integer region = Region;
        algorithm
          state := ThermodynamicState(d = density_pT(p, T, region = region), T = T, phase = 1, h = specificEnthalpy_pT(p, T, region = region), p = p);
          annotation(Inline = true);
        end setState_pTX;
      end WaterIF97_base;

      package IF97_Utilities
        package BaseIF97
          record IterationData
            constant Integer IMAX = 50;
            constant Real DELP = 1.0e-6;
            constant Real DELS = 1.0e-8;
            constant Real DELH = 1.0e-8;
            constant Real DELD = 1.0e-8;
          end IterationData;

          record data
            constant SI.SpecificHeatCapacity RH2O = 461.526;
            constant SI.MolarMass MH2O = 0.01801528;
            constant SI.Temperature TSTAR1 = 1386.0;
            constant SI.Pressure PSTAR1 = 16.53e6;
            constant SI.Temperature TSTAR2 = 540.0;
            constant SI.Pressure PSTAR2 = 1.0e6;
            constant SI.Temperature TSTAR5 = 1000.0;
            constant SI.Pressure PSTAR5 = 1.0e6;
            constant SI.SpecificEnthalpy HSTAR1 = 2.5e6;
            constant Real IPSTAR = 1.0e-6;
            constant Real IHSTAR = 5.0e-7;
            constant SI.Temperature TLIMIT1 = 623.15;
            constant SI.Temperature TLIMIT2 = 1073.15;
            constant SI.Temperature TLIMIT5 = 2273.15;
            constant SI.Pressure PLIMIT1 = 100.0e6;
            constant SI.Pressure PLIMIT4A = 16.5292e6;
            constant SI.Pressure PLIMIT5 = 10.0e6;
            constant SI.Pressure PCRIT = 22064000.0;
            constant SI.Temperature TCRIT = 647.096;
            constant SI.Density DCRIT = 322.0;
            constant SI.SpecificEntropy SCRIT = 4412.02148223476;
            constant SI.SpecificEnthalpy HCRIT = 2087546.84511715;
            constant Real[5] n = array(0.34805185628969e3, -0.11671859879975e1, 0.10192970039326e-2, 0.57254459862746e3, 0.13918839778870e2);
          end data;

          record triple
            constant SI.Temperature Ttriple = 273.16;
            constant SI.Pressure ptriple = 611.657;
            constant SI.Density dltriple = 999.792520031617642;
            constant SI.Density dvtriple = 0.485457572477861372e-2;
          end triple;

          package Regions
            function boundary23ofT
              input SI.Temperature t;
              output SI.Pressure p;
            protected
              constant Real[5] n = data.n;
            algorithm
              p := 1.0e6*(n[1] + t*(n[2] + t*n[3]));
            end boundary23ofT;

            function boundary23ofp
              input SI.Pressure p;
              output SI.Temperature t;
            protected
              constant Real[5] n = data.n;
              Real pi;
            algorithm
              pi := p/1.0e6;
              assert(p > triple.ptriple, "IF97 medium function boundary23ofp called with too low pressure\n" + "p = " + String(p) + " Pa <= " + String(triple.ptriple) + " Pa (triple point pressure)");
              t := n[4] + ((pi - n[5])/n[3])^0.5;
            end boundary23ofp;

            function hlowerofp5
              input SI.Pressure p;
              output SI.SpecificEnthalpy h;
            protected
              Real pi;
            algorithm
              pi := p/data.PSTAR5;
              assert(p > triple.ptriple, "IF97 medium function hlowerofp5 called with too low pressure\n" + "p = " + String(p) + " Pa <= " + String(triple.ptriple) + " Pa (triple point pressure)");
              h := 461526.*(9.01505286876203 + pi*(-0.00979043490246092 + (-0.0000203245575263501 + 3.36540214679088e-7*pi)*pi));
            end hlowerofp5;

            function hupperofp5
              input SI.Pressure p;
              output SI.SpecificEnthalpy h;
            protected
              Real pi;
            algorithm
              pi := p/data.PSTAR5;
              assert(p > triple.ptriple, "IF97 medium function hupperofp5 called with too low pressure\n" + "p = " + String(p) + " Pa <= " + String(triple.ptriple) + " Pa (triple point pressure)");
              h := 461526.*(15.9838891400332 + pi*(-0.000489898813722568 + (-5.01510211858761e-8 + 7.5006972718273e-8*pi)*pi));
            end hupperofp5;

            function hlowerofp1
              input SI.Pressure p;
              output SI.SpecificEnthalpy h;
            protected
              Real pi1;
              Real[3] o;
            algorithm
              pi1 := 7.1 - p/data.PSTAR1;
              assert(p > triple.ptriple, "IF97 medium function hlowerofp1 called with too low pressure\n" + "p = " + String(p) + " Pa <= " + String(triple.ptriple) + " Pa (triple point pressure)");
              o[1] := pi1*pi1;
              o[2] := o[1]*o[1];
              o[3] := o[2]*o[2];
              h := 639675.036*(0.173379420894777 + pi1*(-0.022914084306349 + pi1*(-0.00017146768241932 + pi1*(-4.18695814670391e-6 + pi1*(-2.41630417490008e-7 + pi1*(1.73545618580828e-11 + o[1]*pi1*(8.43755552264362e-14 + o[2]*o[3]*pi1*(5.35429206228374e-35 + o[1]*(-8.12140581014818e-38 + o[1]*o[2]*(-1.43870236842915e-44 + pi1*(1.73894459122923e-45 + (-7.06381628462585e-47 + 9.64504638626269e-49*pi1)*pi1)))))))))));
            end hlowerofp1;

            function hupperofp1
              input SI.Pressure p;
              output SI.SpecificEnthalpy h;
            protected
              Real pi1;
              Real[3] o;
            algorithm
              pi1 := 7.1 - p/data.PSTAR1;
              assert(p > triple.ptriple, "IF97 medium function hupperofp1 called with too low pressure\n" + "p = " + String(p) + " Pa <= " + String(triple.ptriple) + " Pa (triple point pressure)");
              o[1] := pi1*pi1;
              o[2] := o[1]*o[1];
              o[3] := o[2]*o[2];
              h := 639675.036*(2.42896927729349 + pi1*(-0.00141131225285294 + pi1*(0.00143759406818289 + pi1*(0.000125338925082983 + pi1*(0.0000123617764767172 + pi1*(3.17834967400818e-6 + o[1]*pi1*(1.46754947271665e-8 + o[2]*o[3]*pi1*(1.86779322717506e-17 + o[1]*(-4.18568363667416e-19 + o[1]*o[2]*(-9.19148577641497e-22 + pi1*(4.27026404402408e-22 + (-6.66749357417962e-23 + 3.49930466305574e-24*pi1)*pi1)))))))))));
            end hupperofp1;

            function hlowerofp2
              input SI.Pressure p;
              output SI.SpecificEnthalpy h;
            protected
              Real pi;
              Real q1;
              Real q2;
              Real[18] o;
            algorithm
              pi := p/data.PSTAR2;
              assert(p > triple.ptriple, "IF97 medium function hlowerofp2 called with too low pressure\n" + "p = " + String(p) + " Pa <= " + String(triple.ptriple) + " Pa (triple point pressure)");
              q1 := 572.54459862746 + 31.3220101646784*(-13.91883977887 + pi)^0.5;
              q2 := -0.5 + 540./q1;
              o[1] := q1*q1;
              o[2] := o[1]*o[1];
              o[3] := o[2]*o[2];
              o[4] := pi*pi;
              o[5] := o[4]*o[4];
              o[6] := q2*q2;
              o[7] := o[6]*o[6];
              o[8] := o[6]*o[7];
              o[9] := o[5]*o[5];
              o[10] := o[7]*o[7];
              o[11] := o[9]*o[9];
              o[12] := o[10]*o[10];
              o[13] := o[12]*o[12];
              o[14] := o[7]*q2;
              o[15] := o[6]*q2;
              o[16] := o[10]*o[6];
              o[17] := o[13]*o[6];
              o[18] := o[13]*o[6]*q2;
              h := (4.63697573303507e9 + 3.74686560065793*o[2] + 3.57966647812489e-6*o[1]*o[2] + 2.81881548488163e-13*o[3] - 7.64652332452145e7*q1 - 0.00450789338787835*o[2]*q1 - 1.55131504410292e-9*o[1]*o[2]*q1 + o[1]*(2.51383707870341e6 - 4.78198198764471e6*o[10]*o[11]*o[12]*o[13]*o[4] + 49.9651389369988*o[11]*o[12]*o[13]*o[4]*o[5]*o[7] + o[15]*o[4]*(1.03746636552761e-13 - 0.00349547959376899*o[16] - 2.55074501962569e-7*o[8])*o[9] + (-242662.235426958*o[10]*o[12] - 3.46022402653609*o[16])*o[4]*o[5]*pi + o[4]*(0.109336249381227 - 2248.08924686956*o[14] - 354742.725841972*o[17] - 24.1331193696374*o[6])*pi - 3.09081828396912e-19*o[11]*o[12]*o[5]*o[7]*pi - 1.24107527851371e-8*o[11]*o[13]*o[4]*o[5]*o[6]*o[7]*pi + 3.99891272904219*o[5]*o[8]*pi + 0.0641817365250892*o[10]*o[7]*o[9]*pi + pi*(-4444.87643334512 - 75253.6156722047*o[14] - 43051.9020511789*o[6] - 22926.6247146068*q2) + o[4]*(-8.23252840892034 - 3927.0508365636*o[15] - 239.325789467604*o[18] - 76407.3727417716*o[8] - 94.4508644545118*q2) + 0.360567666582363*o[5]*(-0.0161221195808321 + q2)*(0.0338039844460968 + q2) + o[11]*(-0.000584580992538624*o[10]*o[12]*o[7] + 1.33248030241755e6*o[12]*o[13]*q2) + o[9]*(-7.38502736990986e7*o[18] + 0.0000224425477627799*o[6]*o[7]*q2) + o[4]*o[5]*(-2.08438767026518e8*o[17] - 0.0000124971648677697*o[6] - 8442.30378348203*o[10]*o[6]*o[7]*q2) + o[11]*o[9]*(4.73594929247646e-22*o[10]*o[12]*q2 - 13.6411358215175*o[10]*o[12]*o[13]*q2 + 5.52427169406836e-10*o[13]*o[6]*o[7]*q2) + o[11]*o[5]*(2.67174673301715e-6*o[17] + 4.44545133805865e-18*o[12]*o[6]*q2 - 50.2465185106411*o[10]*o[13]*o[6]*o[7]*q2)))/o[1];
            end hlowerofp2;

            function hupperofp2
              input SI.Pressure p;
              output SI.SpecificEnthalpy h;
            protected
              Real pi;
              Real[2] o;
            algorithm
              pi := p/data.PSTAR2;
              assert(p > triple.ptriple, "IF97 medium function hupperofp2 called with too low pressure\n" + "p = " + String(p) + " Pa <= " + String(triple.ptriple) + " Pa (triple point pressure)");
              o[1] := pi*pi;
              o[2] := o[1]*o[1]*o[1];
              h := 4.16066337647071e6 + pi*(-4518.48617188327 + pi*(-8.53409968320258 + pi*(0.109090430596056 + pi*(-0.000172486052272327 + pi*(4.2261295097284e-15 + pi*(-1.27295130636232e-10 + pi*(-3.79407294691742e-25 + pi*(7.56960433802525e-23 + pi*(7.16825117265975e-32 + pi*(3.37267475986401e-21 + (-7.5656940729795e-74 + o[1]*(-8.00969737237617e-134 + (1.6746290980312e-65 + pi*(-3.71600586812966e-69 + pi*(8.06630589170884e-129 + (-1.76117969553159e-103 + 1.88543121025106e-84*pi)*pi)))*o[1]))*o[2]))))))))));
            end hupperofp2;

            function d1n
              input SI.Pressure p;
              input SI.Temperature T;
              output SI.Density d;
            protected
              Real pi;
              Real pi1;
              Real tau;
              Real tau1;
              Real gpi;
              Real[11] o;
            algorithm
              pi := p/data.PSTAR1;
              tau := data.TSTAR1/T;
              pi1 := 7.1 - pi;
              tau1 := tau - 1.222;
              o[1] := tau1*tau1;
              o[2] := o[1]*o[1];
              o[3] := o[2]*o[2];
              o[4] := o[1]*o[2];
              o[5] := o[1]*tau1;
              o[6] := o[2]*tau1;
              o[7] := pi1*pi1;
              o[8] := o[7]*o[7];
              o[9] := o[8]*o[8];
              o[10] := o[3]*o[3];
              o[11] := o[10]*o[10];
              gpi := pi1*(pi1*((0.000095038934535162 + o[2]*(8.4812393955936e-6 + 2.55615384360309e-9*o[4]))/o[2] + pi1*((8.9701127632e-6 + (2.60684891582404e-6 + 5.7366919751696e-13*o[2]*o[3])*o[5])/o[6] + pi1*(2.02584984300585e-6/o[3] + o[7]*pi1*(o[8]*o[9]*pi1*(o[7]*(o[7]*o[8]*(-7.63737668221055e-22/(o[1]*o[11]*o[2]) + pi1*(pi1*(-5.65070932023524e-23/(o[11]*o[3]) + (2.99318679335866e-24*pi1)/(o[11]*o[3]*tau1)) + 3.5842867920213e-22/(o[1]*o[11]*o[2]*tau1))) - 3.33001080055983e-19/(o[1]*o[10]*o[2]*o[3]*tau1)) + 1.44400475720615e-17/(o[10]*o[2]*o[3]*tau1)) + (1.01874413933128e-8 + 1.39398969845072e-9*o[6])/(o[1]*o[3]*tau1))))) + (0.00094368642146534 + o[5]*(0.00060003561586052 + (-0.000095322787813974 + o[1]*(8.8283690661692e-6 + 1.45389992595188e-15*o[1]*o[2]*o[3]))*tau1))/o[5]) + (-0.00028319080123804 + o[1]*(0.00060706301565874 + o[4]*(0.018990068218419 + tau1*(0.032529748770505 + (0.021841717175414 + 0.00005283835796993*o[1])*tau1))))/(o[3]*tau1);
              d := p/(data.RH2O*T*pi*gpi);
            end d1n;

            function d2n
              input SI.Pressure p;
              input SI.Temperature T;
              output SI.Density d;
            protected
              Real pi;
              Real tau;
              Real tau2;
              Real gpi;
              Real[12] o;
            algorithm
              pi := p/data.PSTAR2;
              tau := data.TSTAR2/T;
              tau2 := tau - 0.5;
              o[1] := tau2*tau2;
              o[2] := o[1]*tau2;
              o[3] := o[1]*o[1];
              o[4] := o[3]*o[3];
              o[5] := o[4]*o[4];
              o[6] := o[3]*o[4]*o[5]*tau2;
              o[7] := o[3]*o[4]*tau2;
              o[8] := o[1]*o[3]*o[4];
              o[9] := pi*pi;
              o[10] := o[9]*o[9];
              o[11] := o[3]*o[5]*tau2;
              o[12] := o[5]*o[5];
              gpi := (1. + pi*(-0.0017731742473213 + tau2*(-0.017834862292358 + tau2*(-0.045996013696365 + (-0.057581259083432 - 0.05032527872793*o[2])*tau2)) + pi*(tau2*(-0.000066065283340406 + (-0.0003789797503263 + o[1]*(-0.007878555448671 + o[2]*(-0.087594591301146 - 0.000053349095828174*o[6])))*tau2) + pi*(6.1445213076927e-8 + (1.31612001853305e-6 + o[1]*(-0.00009683303171571 + o[2]*(-0.0045101773626444 - 0.122004760687947*o[6])))*tau2 + pi*(tau2*(-3.15389238237468e-9 + (5.116287140914e-8 + 1.92901490874028e-6*tau2)*tau2) + pi*(0.0000114610381688305*o[1]*o[3]*tau2 + pi*(o[2]*(-1.00288598706366e-10 + o[7]*(-0.012702883392813 - 143.374451604624*o[1]*o[5]*tau2)) + pi*(-4.1341695026989e-17 + o[1]*o[4]*(-8.8352662293707e-6 - 0.272627897050173*o[8])*tau2 + pi*(o[4]*(9.0049690883672e-11 - 65.8490727183984*o[3]*o[4]*o[5]) + pi*(1.78287415218792e-7*o[7] + pi*(o[3]*(1.0406965210174e-18 + o[1]*(-1.0234747095929e-12 - 1.0018179379511e-8*o[3])*o[3]) + o[10]*o[9]*((-1.29412653835176e-9 + 1.71088510070544*o[11])*o[6] + o[9]*(-6.05920510335078*o[12]*o[4]*o[5]*tau2 + o[9]*(o[3]*o[5]*(1.78371690710842e-23 + o[1]*o[3]*o[4]*(6.1258633752464e-12 - 0.000084004935396416*o[7])*tau2) + pi*(-1.24017662339842e-24*o[11] + pi*(0.0000832192847496054*o[12]*o[3]*o[5]*tau2 + pi*(o[1]*o[4]*o[5]*(1.75410265428146e-27 + (1.32995316841867e-15 - 0.0000226487297378904*o[1]*o[5])*o[8])*pi - 2.93678005497663e-14*o[1]*o[12]*o[3]*tau2)))))))))))))))))/pi;
              d := p/(data.RH2O*T*pi*gpi);
            end d2n;

            function hl_p_R4b
              input SI.Pressure p;
              output SI.SpecificEnthalpy h;
            protected
              Real x;
            algorithm
              x := Modelica.Math.acos(p/data.PCRIT);
              h := (1 + x*(-0.4945586958175176 + x*(1.346800016564904 + x*(-3.889388153209752 + x*(6.679385472887931 + x*(-6.75820241066552 + x*(3.558919744656498 + (-0.7179818554978939 - 0.0001152032945617821*x)*x)))))))*data.HCRIT;
            end hl_p_R4b;

            function hv_p_R4b
              input SI.Pressure p;
              output SI.SpecificEnthalpy h;
            protected
              Real x;
            algorithm
              x := Modelica.Math.acos(p/data.PCRIT);
              h := (1 + x*(0.4880153718655694 + x*(0.2079670746250689 + x*(-6.084122698421623 + x*(25.08887602293532 + x*(-48.38215180269516 + x*(45.66489164833212 + (-16.98555442961553 + 0.0006616936460057691*x)*x)))))))*data.HCRIT;
            end hv_p_R4b;

            function rhol_p_R4b
              input SI.Pressure p;
              output SI.Density dl;
            protected
              Real x;
            algorithm
              if (p < data.PCRIT) then
                x := Modelica.Math.acos(p/data.PCRIT);
                dl := (1 + x*(1.903224079094824 + x*(-2.5314861802401123 + x*(-8.191449323843552 + x*(94.34196116778385 + x*(-369.3676833623383 + x*(796.6627910598293 + x*(-994.5385383600702 + x*(673.2581177021598 + (-191.43077336405156 + 0.00052536560808895*x)*x)))))))))*data.DCRIT;
              else
                dl := data.DCRIT;
              end if;
            end rhol_p_R4b;

            function rhov_p_R4b
              input SI.Pressure p;
              output SI.Density dv;
            protected
              Real x;
            algorithm
              if (p < data.PCRIT) then
                x := Modelica.Math.acos(p/data.PCRIT);
                dv := (1 + x*(-1.8463850803362596 + x*(-1.1447872718878493 + x*(59.18702203076563 + x*(-403.5391431811611 + x*(1437.2007245332388 + x*(-3015.853540307519 + x*(3740.5790348670057 + x*(-2537.375817253895 + (725.8761975803782 - 0.0011151111658332337*x)*x)))))))))*data.DCRIT;
              else
                dv := data.DCRIT;
              end if;
            end rhov_p_R4b;

            function boilingcurve_p
              input SI.Pressure p;
              output Common.IF97PhaseBoundaryProperties bpro;
            protected
              Common.GibbsDerivs g;
              Common.HelmholtzDerivs f;
              SI.Pressure plim = min(p, data.PCRIT - 1e-7);
            algorithm
              bpro.R_s := data.RH2O;
              bpro.T := Basic.tsat(plim);
              bpro.dpT := Basic.dptofT(bpro.T);
              bpro.region3boundary := bpro.T > data.TLIMIT1;
              if not bpro.region3boundary then
                g := Basic.g1(p, bpro.T);
                bpro.d := p/(bpro.R_s*bpro.T*g.pi*g.gpi);
                bpro.h := if p > plim then data.HCRIT else bpro.R_s*bpro.T*g.tau*g.gtau;
                bpro.s := g.R_s*(g.tau*g.gtau - g.g);
                bpro.cp := -bpro.R_s*g.tau*g.tau*g.gtautau;
                bpro.vt := bpro.R_s/p*(g.pi*g.gpi - g.tau*g.pi*g.gtaupi);
                bpro.vp := bpro.R_s*bpro.T/(p*p)*g.pi*g.pi*g.gpipi;
                bpro.pt := -p/bpro.T*(g.gpi - g.tau*g.gtaupi)/(g.gpipi*g.pi);
                bpro.pd := -bpro.R_s*bpro.T*g.gpi*g.gpi/(g.gpipi);
              else
                bpro.d := rhol_p_R4b(plim);
                f := Basic.f3(bpro.d, bpro.T);
                bpro.h := hl_p_R4b(plim);
                bpro.s := f.R_s*(f.tau*f.ftau - f.f);
                bpro.cv := bpro.R_s*(-f.tau*f.tau*f.ftautau);
                bpro.pt := bpro.R_s*bpro.d*f.delta*(f.fdelta - f.tau*f.fdeltatau);
                bpro.pd := bpro.R_s*bpro.T*f.delta*(2.0*f.fdelta + f.delta*f.fdeltadelta);
              end if;
            end boilingcurve_p;

            function dewcurve_p
              input SI.Pressure p;
              output Common.IF97PhaseBoundaryProperties bpro;
            protected
              Common.GibbsDerivs g;
              Common.HelmholtzDerivs f;
              SI.Pressure plim = min(p, data.PCRIT - 1e-7);
            algorithm
              bpro.R_s := data.RH2O;
              bpro.T := Basic.tsat(plim);
              bpro.dpT := Basic.dptofT(bpro.T);
              bpro.region3boundary := bpro.T > data.TLIMIT1;
              if not bpro.region3boundary then
                g := Basic.g2(p, bpro.T);
                bpro.d := p/(bpro.R_s*bpro.T*g.pi*g.gpi);
                bpro.h := if p > plim then data.HCRIT else bpro.R_s*bpro.T*g.tau*g.gtau;
                bpro.s := g.R_s*(g.tau*g.gtau - g.g);
                bpro.cp := -bpro.R_s*g.tau*g.tau*g.gtautau;
                bpro.vt := bpro.R_s/p*(g.pi*g.gpi - g.tau*g.pi*g.gtaupi);
                bpro.vp := bpro.R_s*bpro.T/(p*p)*g.pi*g.pi*g.gpipi;
                bpro.pt := -p/bpro.T*(g.gpi - g.tau*g.gtaupi)/(g.gpipi*g.pi);
                bpro.pd := -bpro.R_s*bpro.T*g.gpi*g.gpi/(g.gpipi);
              else
                bpro.d := rhov_p_R4b(plim);
                f := Basic.f3(bpro.d, bpro.T);
                bpro.h := hv_p_R4b(plim);
                bpro.s := f.R_s*(f.tau*f.ftau - f.f);
                bpro.cv := bpro.R_s*(-f.tau*f.tau*f.ftautau);
                bpro.pt := bpro.R_s*bpro.d*f.delta*(f.fdelta - f.tau*f.fdeltatau);
                bpro.pd := bpro.R_s*bpro.T*f.delta*(2.0*f.fdelta + f.delta*f.fdeltadelta);
              end if;
            end dewcurve_p;

            function hvl_p
              input SI.Pressure p;
              input Common.IF97PhaseBoundaryProperties bpro;
              output SI.SpecificEnthalpy h;
            algorithm
              h := bpro.h;
              annotation(derivative(noDerivative = bpro) = hvl_p_der, Inline = false, LateInline = true);
            end hvl_p;

            function hl_p
              input SI.Pressure p;
              output SI.SpecificEnthalpy h;
            algorithm
              h := hvl_p(p, boilingcurve_p(p));
              annotation(Inline = true);
            end hl_p;

            function hv_p
              input SI.Pressure p;
              output SI.SpecificEnthalpy h;
            algorithm
              h := hvl_p(p, dewcurve_p(p));
              annotation(Inline = true);
            end hv_p;

            function hvl_p_der
              input SI.Pressure p;
              input Common.IF97PhaseBoundaryProperties bpro;
              input Real p_der;
              output Real h_der;
            algorithm
              if bpro.region3boundary then
                h_der := ((bpro.d*bpro.pd - bpro.T*bpro.pt)*p_der + (bpro.T*bpro.pt*bpro.pt + bpro.d*bpro.d*bpro.pd*bpro.cv)/bpro.dpT*p_der)/(bpro.pd*bpro.d*bpro.d);
              else
                h_der := (1/bpro.d - bpro.T*bpro.vt)*p_der + bpro.cp/bpro.dpT*p_der;
              end if;
              annotation(Inline = true);
            end hvl_p_der;

            function rhovl_p
              input SI.Pressure p;
              input Common.IF97PhaseBoundaryProperties bpro;
              output SI.Density rho;
            algorithm
              rho := bpro.d;
              annotation(derivative(noDerivative = bpro) = rhovl_p_der, Inline = false, LateInline = true);
            end rhovl_p;

            function rhol_p
              input SI.Pressure p;
              output SI.Density rho;
            algorithm
              rho := rhovl_p(p, boilingcurve_p(p));
              annotation(Inline = true);
            end rhol_p;

            function rhov_p
              input SI.Pressure p;
              output SI.Density rho;
            algorithm
              rho := rhovl_p(p, dewcurve_p(p));
              annotation(Inline = true);
            end rhov_p;

            function rhovl_p_der
              input SI.Pressure p;
              input Common.IF97PhaseBoundaryProperties bpro;
              input Real p_der;
              output Real d_der;
            algorithm
              d_der := if bpro.region3boundary then (p_der - bpro.pt*p_der/bpro.dpT)/bpro.pd else -bpro.d*bpro.d*(bpro.vp + bpro.vt/bpro.dpT)*p_der;
              annotation(Inline = true);
            end rhovl_p_der;

            function rhol_T
              input SI.Temperature T;
              output SI.Density d;
            protected
              SI.Pressure p;
            algorithm
              p := Basic.psat(T);
              if T < data.TLIMIT1 then
                d := d1n(p, T);
              elseif T < data.TCRIT then
                d := rhol_p_R4b(p);
              else
                d := data.DCRIT;
              end if;
            end rhol_T;

            function rhov_T
              input SI.Temperature T;
              output SI.Density d;
            protected
              SI.Pressure p;
            algorithm
              p := Basic.psat(T);
              if T < data.TLIMIT1 then
                d := d2n(p, T);
              elseif T < data.TCRIT then
                d := rhov_p_R4b(p);
              else
                d := data.DCRIT;
              end if;
            end rhov_T;

            function region_ph
              input SI.Pressure p;
              input SI.SpecificEnthalpy h;
              input Integer phase = 0;
              input Integer mode = 0;
              output Integer region;
            protected
              Boolean hsubcrit;
              SI.Temperature Ttest;
              constant Real[5] n = data.n;
              SI.SpecificEnthalpy hl;
              SI.SpecificEnthalpy hv;
            algorithm
              if (mode <> 0) then
                region := mode;
              else
                hl := hl_p(p);
                hv := hv_p(p);
                if (phase == 2) then
                  region := 4;
                else
                  if (p < triple.ptriple) or (p > data.PLIMIT1) or (h < hlowerofp1(p)) or ((p < 10.0e6) and (h > hupperofp5(p))) or ((p >= 10.0e6) and (h > hupperofp2(p))) then
                    region := -1;
                  else
                    hsubcrit := (h < data.HCRIT);
                    if (p < data.PLIMIT4A) then
                      if hsubcrit then
                        if (phase == 1) then
                          region := 1;
                        else
                          if (h < Isentropic.hofpT1(p, Basic.tsat(p))) then
                            region := 1;
                          else
                            region := 4;
                          end if;
                        end if;
                      else
                        if (h > hlowerofp5(p)) then
                          if ((p < data.PLIMIT5) and (h < hupperofp5(p))) then
                            region := 5;
                          else
                            region := -2;
                          end if;
                        else
                          if (phase == 1) then
                            region := 2;
                          else
                            if (h > Isentropic.hofpT2(p, Basic.tsat(p))) then
                              region := 2;
                            else
                              region := 4;
                            end if;
                          end if;
                        end if;
                      end if;
                    else
                      if hsubcrit then
                        if h < hupperofp1(p) then
                          region := 1;
                        else
                          if h < hl or p > data.PCRIT then
                            region := 3;
                          else
                            region := 4;
                          end if;
                        end if;
                      else
                        if (h > hlowerofp2(p)) then
                          region := 2;
                        else
                          if h > hv or p > data.PCRIT then
                            region := 3;
                          else
                            region := 4;
                          end if;
                        end if;
                      end if;
                    end if;
                  end if;
                end if;
              end if;
            end region_ph;

            function region_pT
              input SI.Pressure p;
              input SI.Temperature T;
              input Integer mode = 0;
              output Integer region;
            algorithm
              if (mode <> 0) then
                region := mode;
              else
                if p < data.PLIMIT4A then
                  if T > data.TLIMIT2 then
                    region := 5;
                  elseif T > Basic.tsat(p) then
                    region := 2;
                  else
                    region := 1;
                  end if;
                else
                  if T < data.TLIMIT1 then
                    region := 1;
                  elseif T < boundary23ofp(p) then
                    region := 3;
                  else
                    region := 2;
                  end if;
                end if;
              end if;
            end region_pT;

            function region_dT
              input SI.Density d;
              input SI.Temperature T;
              input Integer phase = 0;
              input Integer mode = 0;
              output Integer region;
            protected
              Boolean Tovercrit;
              SI.Pressure p23;
            algorithm
              Tovercrit := T > data.TCRIT;
              if (mode <> 0) then
                region := mode;
              else
                p23 := boundary23ofT(T);
                if T > data.TLIMIT2 then
                  if d < 20.5655874106483 then
                    region := 5;
                  else
                    assert(false, "Out of valid region for IF97, pressure above region 5!");
                  end if;
                elseif Tovercrit then
                  if d > d2n(p23, T) and T > data.TLIMIT1 then
                    region := 3;
                  elseif T < data.TLIMIT1 then
                    region := 1;
                  else
                    region := 2;
                  end if;
                elseif (d > rhol_T(T)) then
                  if T < data.TLIMIT1 then
                    region := 1;
                  else
                    region := 3;
                  end if;
                elseif (d < rhov_T(T)) then
                  if (d > d2n(p23, T) and T > data.TLIMIT1) then
                    region := 3;
                  else
                    region := 2;
                  end if;
                else
                  region := 4;
                end if;
              end if;
            end region_dT;
          end Regions;

          package Basic
            function g1
              input SI.Pressure p;
              input SI.Temperature T;
              output Modelica.Media.Common.GibbsDerivs g;
            protected
              Real pi1;
              Real tau1;
              Real[45] o;
              Real pl;
            algorithm
              pl := min(p, data.PCRIT - 1);
              assert(p > triple.ptriple, "IF97 medium function g1 called with too low pressure\n" + "p = " + String(p) + " Pa <= " + String(triple.ptriple) + " Pa (triple point pressure)");
              assert(p <= 100.0e6, "IF97 medium function g1: the input pressure (= " + String(p) + " Pa) is higher than 100 MPa");
              assert(T >= 273.15, "IF97 medium function g1: the temperature (= " + String(T) + " K) is lower than 273.15 K!");
              g.p := p;
              g.T := T;
              g.R_s := data.RH2O;
              g.pi := p/data.PSTAR1;
              g.tau := data.TSTAR1/T;
              pi1 := 7.1000000000000 - g.pi;
              tau1 := -1.22200000000000 + g.tau;
              o[1] := tau1*tau1;
              o[2] := o[1]*o[1];
              o[3] := o[2]*o[2];
              o[4] := o[3]*tau1;
              o[5] := 1/o[4];
              o[6] := o[1]*o[2];
              o[7] := o[1]*tau1;
              o[8] := 1/o[7];
              o[9] := o[1]*o[2]*o[3];
              o[10] := 1/o[2];
              o[11] := o[2]*tau1;
              o[12] := 1/o[11];
              o[13] := o[2]*o[3];
              o[14] := 1/o[3];
              o[15] := pi1*pi1;
              o[16] := o[15]*pi1;
              o[17] := o[15]*o[15];
              o[18] := o[17]*o[17];
              o[19] := o[17]*o[18]*pi1;
              o[20] := o[15]*o[17];
              o[21] := o[3]*o[3];
              o[22] := o[21]*o[21];
              o[23] := o[22]*o[3]*tau1;
              o[24] := 1/o[23];
              o[25] := o[22]*o[3];
              o[26] := 1/o[25];
              o[27] := o[1]*o[2]*o[22]*tau1;
              o[28] := 1/o[27];
              o[29] := o[1]*o[2]*o[22];
              o[30] := 1/o[29];
              o[31] := o[1]*o[2]*o[21]*o[3]*tau1;
              o[32] := 1/o[31];
              o[33] := o[2]*o[21]*o[3]*tau1;
              o[34] := 1/o[33];
              o[35] := o[1]*o[3]*tau1;
              o[36] := 1/o[35];
              o[37] := o[1]*o[3];
              o[38] := 1/o[37];
              o[39] := 1/o[6];
              o[40] := o[1]*o[22]*o[3];
              o[41] := 1/o[40];
              o[42] := 1/o[22];
              o[43] := o[1]*o[2]*o[21]*o[3];
              o[44] := 1/o[43];
              o[45] := 1/o[13];
              g.g := pi1*(pi1*(pi1*(o[10]*(-0.000031679644845054 + o[2]*(-2.82707979853120e-6 - 8.5205128120103e-10*o[6])) + pi1*(o[12]*(-2.24252819080000e-6 + (-6.5171222895601e-7 - 1.43417299379240e-13*o[13])*o[7]) + pi1*(-4.0516996860117e-7*o[14] + o[16]*((-1.27343017416410e-9 - 1.74248712306340e-10*o[11])*o[36] + o[19]*(-6.8762131295531e-19*o[34] + o[15]*(1.44783078285210e-20*o[32] + o[20]*(2.63357816627950e-23*o[30] + pi1*(-1.19476226400710e-23*o[28] + pi1*(1.82280945814040e-24*o[26] - 9.3537087292458e-26*o[24]*pi1))))))))) + o[8]*(-0.00047184321073267 + o[7]*(-0.000300017807930260 + (0.000047661393906987 + o[1]*(-4.4141845330846e-6 - 7.2694996297594e-16*o[9]))*tau1))) + o[5]*(0.000283190801238040 + o[1]*(-0.00060706301565874 + o[6]*(-0.0189900682184190 + tau1*(-0.032529748770505 + (-0.0218417171754140 - 0.000052838357969930*o[1])*tau1))))) + (0.146329712131670 + tau1*(-0.84548187169114 + tau1*(-3.7563603672040 + tau1*(3.3855169168385 + tau1*(-0.95791963387872 + tau1*(0.157720385132280 + (-0.0166164171995010 + 0.00081214629983568*tau1)*tau1))))))/o[1];
              g.gpi := pi1*(pi1*(o[10]*(0.000095038934535162 + o[2]*(8.4812393955936e-6 + 2.55615384360309e-9*o[6])) + pi1*(o[12]*(8.9701127632000e-6 + (2.60684891582404e-6 + 5.7366919751696e-13*o[13])*o[7]) + pi1*(2.02584984300585e-6*o[14] + o[16]*((1.01874413933128e-8 + 1.39398969845072e-9*o[11])*o[36] + o[19]*(1.44400475720615e-17*o[34] + o[15]*(-3.3300108005598e-19*o[32] + o[20]*(-7.6373766822106e-22*o[30] + pi1*(3.5842867920213e-22*o[28] + pi1*(-5.6507093202352e-23*o[26] + 2.99318679335866e-24*o[24]*pi1))))))))) + o[8]*(0.00094368642146534 + o[7]*(0.00060003561586052 + (-0.000095322787813974 + o[1]*(8.8283690661692e-6 + 1.45389992595188e-15*o[9]))*tau1))) + o[5]*(-0.000283190801238040 + o[1]*(0.00060706301565874 + o[6]*(0.0189900682184190 + tau1*(0.032529748770505 + (0.0218417171754140 + 0.000052838357969930*o[1])*tau1))));
              g.gpipi := pi1*(o[10]*(-0.000190077869070324 + o[2]*(-0.0000169624787911872 - 5.1123076872062e-9*o[6])) + pi1*(o[12]*(-0.0000269103382896000 + (-7.8205467474721e-6 - 1.72100759255088e-12*o[13])*o[7]) + pi1*(-8.1033993720234e-6*o[14] + o[16]*((-7.1312089753190e-8 - 9.7579278891550e-9*o[11])*o[36] + o[19]*(-2.88800951441230e-16*o[34] + o[15]*(7.3260237612316e-18*o[32] + o[20]*(2.13846547101895e-20*o[30] + pi1*(-1.03944316968618e-20*o[28] + pi1*(1.69521279607057e-21*o[26] - 9.2788790594118e-23*o[24]*pi1))))))))) + o[8]*(-0.00094368642146534 + o[7]*(-0.00060003561586052 + (0.000095322787813974 + o[1]*(-8.8283690661692e-6 - 1.45389992595188e-15*o[9]))*tau1));
              g.gtau := pi1*(o[38]*(-0.00254871721114236 + o[1]*(0.0042494411096112 + (0.0189900682184190 + (-0.0218417171754140 - 0.000158515073909790*o[1])*o[1])*o[6])) + pi1*(o[10]*(0.00141552963219801 + o[2]*(0.000047661393906987 + o[1]*(-0.0000132425535992538 - 1.23581493705910e-14*o[9]))) + pi1*(o[12]*(0.000126718579380216 - 5.1123076872062e-9*o[37]) + pi1*(o[39]*(0.0000112126409540000 + (1.30342445791202e-6 - 1.43417299379240e-12*o[13])*o[7]) + pi1*(3.2413597488094e-6*o[5] + o[16]*((1.40077319158051e-8 + 1.04549227383804e-9*o[11])*o[45] + o[19]*(1.99410180757040e-17*o[44] + o[15]*(-4.4882754268415e-19*o[42] + o[20]*(-1.00075970318621e-21*o[28] + pi1*(4.6595728296277e-22*o[26] + pi1*(-7.2912378325616e-23*o[24] + 3.8350205789908e-24*o[41]*pi1))))))))))) + o[8]*(-0.292659424263340 + tau1*(0.84548187169114 + o[1]*(3.3855169168385 + tau1*(-1.91583926775744 + tau1*(0.47316115539684 + (-0.066465668798004 + 0.0040607314991784*tau1)*tau1)))));
              g.gtautau := pi1*(o[36]*(0.0254871721114236 + o[1]*(-0.033995528876889 + (-0.037980136436838 - 0.00031703014781958*o[2])*o[6])) + pi1*(o[12]*(-0.0056621185287920 + o[6]*(-0.0000264851071985076 - 1.97730389929456e-13*o[9])) + pi1*((-0.00063359289690108 - 2.55615384360309e-8*o[37])*o[39] + pi1*(pi1*(-0.0000291722377392842*o[38] + o[16]*(o[19]*(-5.9823054227112e-16*o[32] + o[15]*(o[20]*(3.9029628424262e-20*o[26] + pi1*(-1.86382913185108e-20*o[24] + pi1*(2.98940751135026e-21*o[41] - (1.61070864317613e-22*pi1)/(o[1]*o[22]*o[3]*tau1)))) + 1.43624813658928e-17/(o[22]*tau1))) + (-1.68092782989661e-7 - 7.3184459168663e-9*o[11])/(o[2]*o[3]*tau1))) + (-0.000067275845724000 + (-3.9102733737361e-6 - 1.29075569441316e-11*o[13])*o[7])/(o[1]*o[2]*tau1))))) + o[10]*(0.87797827279002 + tau1*(-1.69096374338228 + o[7]*(-1.91583926775744 + tau1*(0.94632231079368 + (-0.199397006394012 + 0.0162429259967136*tau1)*tau1))));
              g.gtaupi := o[38]*(0.00254871721114236 + o[1]*(-0.0042494411096112 + (-0.0189900682184190 + (0.0218417171754140 + 0.000158515073909790*o[1])*o[1])*o[6])) + pi1*(o[10]*(-0.00283105926439602 + o[2]*(-0.000095322787813974 + o[1]*(0.0000264851071985076 + 2.47162987411820e-14*o[9]))) + pi1*(o[12]*(-0.00038015573814065 + 1.53369230616185e-8*o[37]) + pi1*(o[39]*(-0.000044850563816000 + (-5.2136978316481e-6 + 5.7366919751696e-12*o[13])*o[7]) + pi1*(-0.0000162067987440468*o[5] + o[16]*((-1.12061855326441e-7 - 8.3639381907043e-9*o[11])*o[45] + o[19]*(-4.1876137958978e-16*o[44] + o[15]*(1.03230334817355e-17*o[42] + o[20]*(2.90220313924001e-20*o[28] + pi1*(-1.39787184888831e-20*o[26] + pi1*(2.26028372809410e-21*o[24] - 1.22720658527705e-22*o[41]*pi1))))))))));
            end g1;

            function g2
              input SI.Pressure p;
              input SI.Temperature T;
              input Boolean checkLimits = true;
              output Modelica.Media.Common.GibbsDerivs g;
            protected
              Real tau2;
              Real[55] o;
            algorithm
              g.p := p;
              g.T := T;
              g.R_s := data.RH2O;
              if checkLimits then
                assert(p > 0.0, "IF97 medium function g2 called with too low pressure\n" + "p = " + String(p) + " Pa <= 0.0 Pa");
                assert(p <= 100.0e6, "IF97 medium function g2: the input pressure (= " + String(p) + " Pa) is higher than 100 MPa");
                assert(T >= 273.15, "IF97 medium function g2: the temperature (= " + String(T) + " K) is lower than 273.15 K!");
                assert(T <= 1073.15, "IF97 medium function g2: the input temperature (= " + String(T) + " K) is higher than the limit of 1073.15 K");
              else
              end if;
              g.pi := p/data.PSTAR2;
              g.tau := data.TSTAR2/T;
              tau2 := -0.5 + g.tau;
              o[1] := tau2*tau2;
              o[2] := o[1]*tau2;
              o[3] := -0.050325278727930*o[2];
              o[4] := -0.057581259083432 + o[3];
              o[5] := o[4]*tau2;
              o[6] := -0.045996013696365 + o[5];
              o[7] := o[6]*tau2;
              o[8] := -0.0178348622923580 + o[7];
              o[9] := o[8]*tau2;
              o[10] := o[1]*o[1];
              o[11] := o[10]*o[10];
              o[12] := o[11]*o[11];
              o[13] := o[10]*o[11]*o[12]*tau2;
              o[14] := o[1]*o[10]*tau2;
              o[15] := o[10]*o[11]*tau2;
              o[16] := o[1]*o[12]*tau2;
              o[17] := o[1]*o[11]*tau2;
              o[18] := o[1]*o[10]*o[11];
              o[19] := o[10]*o[11]*o[12];
              o[20] := o[1]*o[10];
              o[21] := g.pi*g.pi;
              o[22] := o[21]*o[21];
              o[23] := o[21]*o[22];
              o[24] := o[10]*o[12]*tau2;
              o[25] := o[12]*o[12];
              o[26] := o[11]*o[12]*o[25]*tau2;
              o[27] := o[10]*o[12];
              o[28] := o[1]*o[10]*o[11]*tau2;
              o[29] := o[10]*o[12]*o[25]*tau2;
              o[30] := o[1]*o[10]*o[25]*tau2;
              o[31] := o[1]*o[11]*o[12];
              o[32] := o[1]*o[12];
              o[33] := g.tau*g.tau;
              o[34] := o[33]*o[33];
              o[35] := -0.000053349095828174*o[13];
              o[36] := -0.087594591301146 + o[35];
              o[37] := o[2]*o[36];
              o[38] := -0.0078785554486710 + o[37];
              o[39] := o[1]*o[38];
              o[40] := -0.00037897975032630 + o[39];
              o[41] := o[40]*tau2;
              o[42] := -0.000066065283340406 + o[41];
              o[43] := o[42]*tau2;
              o[44] := 5.7870447262208e-6*tau2;
              o[45] := -0.301951672367580*o[2];
              o[46] := -0.172743777250296 + o[45];
              o[47] := o[46]*tau2;
              o[48] := -0.091992027392730 + o[47];
              o[49] := o[48]*tau2;
              o[50] := o[1]*o[11];
              o[51] := o[10]*o[11];
              o[52] := o[11]*o[12]*o[25];
              o[53] := o[10]*o[12]*o[25];
              o[54] := o[1]*o[10]*o[25];
              o[55] := o[11]*o[12]*tau2;
              g.g := g.pi*(-0.00177317424732130 + o[9] + g.pi*(tau2*(-0.000033032641670203 + (-0.000189489875163150 + o[1]*(-0.0039392777243355 + (-0.043797295650573 - 0.0000266745479140870*o[13])*o[2]))*tau2) + g.pi*(2.04817376923090e-8 + (4.3870667284435e-7 + o[1]*(-0.000032277677238570 + (-0.00150339245421480 - 0.040668253562649*o[13])*o[2]))*tau2 + g.pi*(g.pi*(2.29220763376610e-6*o[14] + g.pi*((-1.67147664510610e-11 + o[15]*(-0.00211714723213550 - 23.8957419341040*o[16]))*o[2] + g.pi*(-5.9059564324270e-18 + o[17]*(-1.26218088991010e-6 - 0.038946842435739*o[18]) + g.pi*(o[11]*(1.12562113604590e-11 - 8.2311340897998*o[19]) + g.pi*(1.98097128020880e-8*o[15] + g.pi*(o[10]*(1.04069652101740e-19 + (-1.02347470959290e-13 - 1.00181793795110e-9*o[10])*o[20]) + o[23]*(o[13]*(-8.0882908646985e-11 + 0.106930318794090*o[24]) + o[21]*(-0.33662250574171*o[26] + o[21]*(o[27]*(8.9185845355421e-25 + (3.06293168762320e-13 - 4.2002467698208e-6*o[15])*o[28]) + g.pi*(-5.9056029685639e-26*o[24] + g.pi*(3.7826947613457e-6*o[29] + g.pi*(-1.27686089346810e-15*o[30] + o[31]*(7.3087610595061e-29 + o[18]*(5.5414715350778e-17 - 9.4369707241210e-7*o[32]))*g.pi)))))))))))) + tau2*(-7.8847309559367e-10 + (1.27907178522850e-8 + 4.8225372718507e-7*tau2)*tau2))))) + (-0.0056087911830200 + g.tau*(0.071452738814550 + g.tau*(-0.40710498239280 + g.tau*(1.42408197144400 + g.tau*(-4.3839511194500 + g.tau*(-9.6927686002170 + g.tau*(10.0866556801800 + (-0.284086326077200 + 0.0212684635330700*g.tau)*g.tau) + Modelica.Math.log(g.pi)))))))/(o[34]*g.tau);
              g.gpi := (1.00000000000000 + g.pi*(-0.00177317424732130 + o[9] + g.pi*(o[43] + g.pi*(6.1445213076927e-8 + (1.31612001853305e-6 + o[1]*(-0.000096833031715710 + (-0.0045101773626444 - 0.122004760687947*o[13])*o[2]))*tau2 + g.pi*(g.pi*(0.0000114610381688305*o[14] + g.pi*((-1.00288598706366e-10 + o[15]*(-0.0127028833928130 - 143.374451604624*o[16]))*o[2] + g.pi*(-4.1341695026989e-17 + o[17]*(-8.8352662293707e-6 - 0.272627897050173*o[18]) + g.pi*(o[11]*(9.0049690883672e-11 - 65.849072718398*o[19]) + g.pi*(1.78287415218792e-7*o[15] + g.pi*(o[10]*(1.04069652101740e-18 + (-1.02347470959290e-12 - 1.00181793795110e-8*o[10])*o[20]) + o[23]*(o[13]*(-1.29412653835176e-9 + 1.71088510070544*o[24]) + o[21]*(-6.0592051033508*o[26] + o[21]*(o[27]*(1.78371690710842e-23 + (6.1258633752464e-12 - 0.000084004935396416*o[15])*o[28]) + g.pi*(-1.24017662339842e-24*o[24] + g.pi*(0.000083219284749605*o[29] + g.pi*(-2.93678005497663e-14*o[30] + o[31]*(1.75410265428146e-27 + o[18]*(1.32995316841867e-15 - 0.0000226487297378904*o[32]))*g.pi)))))))))))) + tau2*(-3.15389238237468e-9 + (5.1162871409140e-8 + 1.92901490874028e-6*tau2)*tau2))))))/g.pi;
              g.gpipi := (-1.00000000000000 + o[21]*(o[43] + g.pi*(1.22890426153854e-7 + (2.63224003706610e-6 + o[1]*(-0.000193666063431420 + (-0.0090203547252888 - 0.244009521375894*o[13])*o[2]))*tau2 + g.pi*(g.pi*(0.000045844152675322*o[14] + g.pi*((-5.0144299353183e-10 + o[15]*(-0.063514416964065 - 716.87225802312*o[16]))*o[2] + g.pi*(-2.48050170161934e-16 + o[17]*(-0.000053011597376224 - 1.63576738230104*o[18]) + g.pi*(o[11]*(6.3034783618570e-10 - 460.94350902879*o[19]) + g.pi*(1.42629932175034e-6*o[15] + g.pi*(o[10]*(9.3662686891566e-18 + (-9.2112723863361e-12 - 9.0163614415599e-8*o[10])*o[20]) + o[23]*(o[13]*(-1.94118980752764e-8 + 25.6632765105816*o[24]) + o[21]*(-103.006486756963*o[26] + o[21]*(o[27]*(3.3890621235060e-22 + (1.16391404129682e-10 - 0.00159609377253190*o[15])*o[28]) + g.pi*(-2.48035324679684e-23*o[24] + g.pi*(0.00174760497974171*o[29] + g.pi*(-6.4609161209486e-13*o[30] + o[31]*(4.0344361048474e-26 + o[18]*(3.05889228736295e-14 - 0.00052092078397148*o[32]))*g.pi)))))))))))) + tau2*(-9.4616771471240e-9 + (1.53488614227420e-7 + o[44])*tau2)))))/o[21];
              g.gtau := (0.0280439559151000 + g.tau*(-0.285810955258200 + g.tau*(1.22131494717840 + g.tau*(-2.84816394288800 + g.tau*(4.3839511194500 + o[33]*(10.0866556801800 + (-0.56817265215440 + 0.063805390599210*g.tau)*g.tau))))))/(o[33]*o[34]) + g.pi*(-0.0178348622923580 + o[49] + g.pi*(-0.000033032641670203 + (-0.00037897975032630 + o[1]*(-0.0157571108973420 + (-0.306581069554011 - 0.00096028372490713*o[13])*o[2]))*tau2 + g.pi*(4.3870667284435e-7 + o[1]*(-0.000096833031715710 + (-0.0090203547252888 - 1.42338887469272*o[13])*o[2]) + g.pi*(-7.8847309559367e-10 + g.pi*(0.0000160454534363627*o[20] + g.pi*(o[1]*(-5.0144299353183e-11 + o[15]*(-0.033874355714168 - 836.35096769364*o[16])) + g.pi*((-0.0000138839897890111 - 0.97367106089347*o[18])*o[50] + g.pi*(o[14]*(9.0049690883672e-11 - 296.320827232793*o[19]) + g.pi*(2.57526266427144e-7*o[51] + g.pi*(o[2]*(4.1627860840696e-19 + (-1.02347470959290e-12 - 1.40254511313154e-8*o[10])*o[20]) + o[23]*(o[19]*(-2.34560435076256e-9 + 5.3465159397045*o[24]) + o[21]*(-19.1874828272775*o[52] + o[21]*(o[16]*(1.78371690710842e-23 + (1.07202609066812e-11 - 0.000201611844951398*o[15])*o[28]) + g.pi*(-1.24017662339842e-24*o[27] + g.pi*(0.000200482822351322*o[53] + g.pi*(-4.9797574845256e-14*o[54] + (1.90027787547159e-27 + o[18]*(2.21658861403112e-15 - 0.000054734430199902*o[32]))*o[55]*g.pi)))))))))))) + (2.55814357045700e-8 + 1.44676118155521e-6*tau2)*tau2))));
              g.gtautau := (-0.168263735490600 + g.tau*(1.42905477629100 + g.tau*(-4.8852597887136 + g.tau*(8.5444918286640 + g.tau*(-8.7679022389000 + o[33]*(-0.56817265215440 + 0.127610781198420*g.tau)*g.tau)))))/(o[33]*o[34]*g.tau) + g.pi*(-0.091992027392730 + (-0.34548755450059 - 1.50975836183790*o[2])*tau2 + g.pi*(-0.00037897975032630 + o[1]*(-0.047271332692026 + (-1.83948641732407 - 0.033609930371750*o[13])*o[2]) + g.pi*((-0.000193666063431420 + (-0.045101773626444 - 48.395221739552*o[13])*o[2])*tau2 + g.pi*(2.55814357045700e-8 + 2.89352236311042e-6*tau2 + g.pi*(0.000096272720618176*o[10]*tau2 + g.pi*((-1.00288598706366e-10 + o[15]*(-0.50811533571252 - 28435.9329015838*o[16]))*tau2 + g.pi*(o[11]*(-0.000138839897890111 - 23.3681054614434*o[18])*tau2 + g.pi*((6.3034783618570e-10 - 10371.2289531477*o[19])*o[20] + g.pi*(3.09031519712573e-6*o[17] + g.pi*(o[1]*(1.24883582522088e-18 + (-9.2112723863361e-12 - 1.82330864707100e-7*o[10])*o[20]) + o[23]*(o[1]*o[11]*o[12]*(-6.5676921821352e-8 + 261.979281045521*o[24])*tau2 + o[21]*(-1074.49903832754*o[1]*o[10]*o[12]*o[25]*tau2 + o[21]*((3.3890621235060e-22 + (3.6448887082716e-10 - 0.0094757567127157*o[15])*o[28])*o[32] + g.pi*(-2.48035324679684e-23*o[16] + g.pi*(0.0104251067622687*o[1]*o[12]*o[25]*tau2 + g.pi*(o[11]*o[12]*(4.7506946886790e-26 + o[18]*(8.6446955947214e-14 - 0.00311986252139440*o[32]))*g.pi - 1.89230784411972e-12*o[10]*o[25]*tau2))))))))))))))));
              g.gtaupi := -0.0178348622923580 + o[49] + g.pi*(-0.000066065283340406 + (-0.00075795950065260 + o[1]*(-0.0315142217946840 + (-0.61316213910802 - 0.00192056744981426*o[13])*o[2]))*tau2 + g.pi*(1.31612001853305e-6 + o[1]*(-0.000290499095147130 + (-0.0270610641758664 - 4.2701666240781*o[13])*o[2]) + g.pi*(-3.15389238237468e-9 + g.pi*(0.000080227267181813*o[20] + g.pi*(o[1]*(-3.00865796119098e-10 + o[15]*(-0.203246134285008 - 5018.1058061618*o[16])) + g.pi*((-0.000097187928523078 - 6.8156974262543*o[18])*o[50] + g.pi*(o[14]*(7.2039752706938e-10 - 2370.56661786234*o[19]) + g.pi*(2.31773639784430e-6*o[51] + g.pi*(o[2]*(4.1627860840696e-18 + (-1.02347470959290e-11 - 1.40254511313154e-7*o[10])*o[20]) + o[23]*(o[19]*(-3.7529669612201e-8 + 85.544255035272*o[24]) + o[21]*(-345.37469089099*o[52] + o[21]*(o[16]*(3.5674338142168e-22 + (2.14405218133624e-10 - 0.0040322368990280*o[15])*o[28]) + g.pi*(-2.60437090913668e-23*o[27] + g.pi*(0.0044106220917291*o[53] + g.pi*(-1.14534422144089e-12*o[54] + (4.5606669011318e-26 + o[18]*(5.3198126736747e-14 - 0.00131362632479764*o[32]))*o[55]*g.pi)))))))))))) + (1.02325742818280e-7 + o[44])*tau2)));
            end g2;

            function f3
              input SI.Density d;
              input SI.Temperature T;
              output Modelica.Media.Common.HelmholtzDerivs f;
            protected
              Real[40] o;
            algorithm
              f.T := T;
              f.d := d;
              f.R_s := data.RH2O;
              f.tau := data.TCRIT/T;
              f.delta := if (d == data.DCRIT and T == data.TCRIT) then 1 - Modelica.Constants.eps else abs(d/data.DCRIT);
              o[1] := f.tau*f.tau;
              o[2] := o[1]*o[1];
              o[3] := o[2]*f.tau;
              o[4] := o[1]*f.tau;
              o[5] := o[2]*o[2];
              o[6] := o[1]*o[5]*f.tau;
              o[7] := o[5]*f.tau;
              o[8] := -0.64207765181607*o[1];
              o[9] := 0.88521043984318 + o[8];
              o[10] := o[7]*o[9];
              o[11] := -1.15244078066810 + o[10];
              o[12] := o[11]*o[2];
              o[13] := -1.26543154777140 + o[12];
              o[14] := o[1]*o[13];
              o[15] := o[1]*o[2]*o[5]*f.tau;
              o[16] := o[2]*o[5];
              o[17] := o[1]*o[5];
              o[18] := o[5]*o[5];
              o[19] := o[1]*o[18]*o[2];
              o[20] := o[1]*o[18]*o[2]*f.tau;
              o[21] := o[18]*o[5];
              o[22] := o[1]*o[18]*o[5];
              o[23] := 0.251168168486160*o[2];
              o[24] := 0.078841073758308 + o[23];
              o[25] := o[15]*o[24];
              o[26] := -6.1005234513930 + o[25];
              o[27] := o[26]*f.tau;
              o[28] := 9.7944563083754 + o[27];
              o[29] := o[2]*o[28];
              o[30] := -1.70429417648412 + o[29];
              o[31] := o[1]*o[30];
              o[32] := f.delta*f.delta;
              o[33] := -10.9153200808732*o[1];
              o[34] := 13.2781565976477 + o[33];
              o[35] := o[34]*o[7];
              o[36] := -6.9146446840086 + o[35];
              o[37] := o[2]*o[36];
              o[38] := -2.53086309554280 + o[37];
              o[39] := o[38]*f.tau;
              o[40] := o[18]*o[5]*f.tau;
              f.f := -15.7328452902390 + f.tau*(20.9443969743070 + (-7.6867707878716 + o[3]*(2.61859477879540 + o[4]*(-2.80807811486200 + o[1]*(1.20533696965170 - 0.0084566812812502*o[6]))))*f.tau) + f.delta*(o[14] + f.delta*(0.38493460186671 + o[1]*(-0.85214708824206 + o[2]*(4.8972281541877 + (-3.05026172569650 + o[15]*(0.039420536879154 + 0.125584084243080*o[2]))*f.tau)) + f.delta*(-0.279993296987100 + o[1]*(1.38997995694600 + o[1]*(-2.01899150235700 + o[16]*(-0.0082147637173963 - 0.47596035734923*o[17]))) + f.delta*(0.043984074473500 + o[1]*(-0.44476435428739 + o[1]*(0.90572070719733 + 0.70522450087967*o[19])) + f.delta*(f.delta*(-0.0221754008730960 + o[1]*(0.094260751665092 + 0.164362784479610*o[21]) + f.delta*(-0.0135033722413480*o[1] + f.delta*(-0.0148343453524720*o[22] + f.delta*(o[1]*(0.00057922953628084 + 0.0032308904703711*o[21]) + f.delta*(0.000080964802996215 - 0.000044923899061815*f.delta*o[22] - 0.000165576797950370*f.tau))))) + (0.107705126263320 + o[1]*(-0.32913623258954 - 0.50871062041158*o[20]))*f.tau))))) + 1.06580700285130*Modelica.Math.log(f.delta);
              f.fdelta := (1.06580700285130 + f.delta*(o[14] + f.delta*(0.76986920373342 + o[31] + f.delta*(-0.83997989096130 + o[1]*(4.1699398708380 + o[1]*(-6.0569745070710 + o[16]*(-0.0246442911521889 - 1.42788107204769*o[17]))) + f.delta*(0.175936297894000 + o[1]*(-1.77905741714956 + o[1]*(3.6228828287893 + 2.82089800351868*o[19])) + f.delta*(f.delta*(-0.133052405238576 + o[1]*(0.56556450999055 + 0.98617670687766*o[21]) + f.delta*(-0.094523605689436*o[1] + f.delta*(-0.118674762819776*o[22] + f.delta*(o[1]*(0.0052130658265276 + 0.0290780142333399*o[21]) + f.delta*(0.00080964802996215 - 0.00049416288967996*f.delta*o[22] - 0.00165576797950370*f.tau))))) + (0.53852563131660 + o[1]*(-1.64568116294770 - 2.54355310205790*o[20]))*f.tau))))))/f.delta;
              f.fdeltadelta := (-1.06580700285130 + o[32]*(0.76986920373342 + o[31] + f.delta*(-1.67995978192260 + o[1]*(8.3398797416760 + o[1]*(-12.1139490141420 + o[16]*(-0.049288582304378 - 2.85576214409538*o[17]))) + f.delta*(0.52780889368200 + o[1]*(-5.3371722514487 + o[1]*(10.8686484863680 + 8.4626940105560*o[19])) + f.delta*(f.delta*(-0.66526202619288 + o[1]*(2.82782254995276 + 4.9308835343883*o[21]) + f.delta*(-0.56714163413662*o[1] + f.delta*(-0.83072333973843*o[22] + f.delta*(o[1]*(0.041704526612220 + 0.232624113866719*o[21]) + f.delta*(0.0072868322696594 - 0.0049416288967996*f.delta*o[22] - 0.0149019118155333*f.tau))))) + (2.15410252526640 + o[1]*(-6.5827246517908 - 10.1742124082316*o[20]))*f.tau)))))/o[32];
              f.ftau := 20.9443969743070 + (-15.3735415757432 + o[3]*(18.3301634515678 + o[4]*(-28.0807811486200 + o[1]*(14.4640436358204 - 0.194503669468755*o[6]))))*f.tau + f.delta*(o[39] + f.delta*(f.tau*(-1.70429417648412 + o[2]*(29.3833689251262 + (-21.3518320798755 + o[15]*(0.86725181134139 + 3.2651861903201*o[2]))*f.tau)) + f.delta*((2.77995991389200 + o[1]*(-8.0759660094280 + o[16]*(-0.131436219478341 - 12.3749692910800*o[17])))*f.tau + f.delta*((-0.88952870857478 + o[1]*(3.6228828287893 + 18.3358370228714*o[19]))*f.tau + f.delta*(0.107705126263320 + o[1]*(-0.98740869776862 - 13.2264761307011*o[20]) + f.delta*((0.188521503330184 + 4.2734323964699*o[21])*f.tau + f.delta*(-0.0270067444826960*f.tau + f.delta*(-0.38569297916427*o[40] + f.delta*(f.delta*(-0.000165576797950370 - 0.00116802137560719*f.delta*o[40]) + (0.00115845907256168 + 0.084003152229649*o[21])*f.tau)))))))));
              f.ftautau := -15.3735415757432 + o[3]*(109.980980709407 + o[4]*(-252.727030337580 + o[1]*(159.104479994024 - 4.2790807283126*o[6]))) + f.delta*(-2.53086309554280 + o[2]*(-34.573223420043 + (185.894192367068 - 174.645121293971*o[1])*o[7]) + f.delta*(-1.70429417648412 + o[2]*(146.916844625631 + (-128.110992479253 + o[15]*(18.2122880381691 + 81.629654758002*o[2]))*f.tau) + f.delta*(2.77995991389200 + o[1]*(-24.2278980282840 + o[16]*(-1.97154329217511 - 309.374232277000*o[17])) + f.delta*(-0.88952870857478 + o[1]*(10.8686484863680 + 458.39592557179*o[19]) + f.delta*(f.delta*(0.188521503330184 + 106.835809911747*o[21] + f.delta*(-0.0270067444826960 + f.delta*(-9.6423244791068*o[21] + f.delta*(0.00115845907256168 + 2.10007880574121*o[21] - 0.0292005343901797*o[21]*o[32])))) + (-1.97481739553724 - 330.66190326753*o[20])*f.tau)))));
              f.fdeltatau := o[39] + f.delta*(f.tau*(-3.4085883529682 + o[2]*(58.766737850252 + (-42.703664159751 + o[15]*(1.73450362268278 + 6.5303723806402*o[2]))*f.tau)) + f.delta*((8.3398797416760 + o[1]*(-24.2278980282840 + o[16]*(-0.39430865843502 - 37.124907873240*o[17])))*f.tau + f.delta*((-3.5581148342991 + o[1]*(14.4915313151573 + 73.343348091486*o[19]))*f.tau + f.delta*(0.53852563131660 + o[1]*(-4.9370434888431 - 66.132380653505*o[20]) + f.delta*((1.13112901998110 + 25.6405943788192*o[21])*f.tau + f.delta*(-0.189047211378872*f.tau + f.delta*(-3.08554383331418*o[40] + f.delta*(f.delta*(-0.00165576797950370 - 0.0128482351316791*f.delta*o[40]) + (0.0104261316530551 + 0.75602837006684*o[21])*f.tau))))))));
            end f3;

            function g5
              input SI.Pressure p;
              input SI.Temperature T;
              output Modelica.Media.Common.GibbsDerivs g;
            protected
              Real[11] o;
            algorithm
              assert(p > 0.0, "IF97 medium function g5 called with too low pressure\n" + "p = " + String(p) + " Pa <=  0.0 Pa");
              assert(p <= data.PLIMIT5, "IF97 medium function g5: input pressure (= " + String(p) + " Pa) is higher than 10 MPa in region 5");
              assert(T <= 2273.15, "IF97 medium function g5: input temperature (= " + String(T) + " K) is higher than limit of 2273.15K in region 5");
              g.p := p;
              g.T := T;
              g.R_s := data.RH2O;
              g.pi := p/data.PSTAR5;
              g.tau := data.TSTAR5/T;
              o[1] := g.tau*g.tau;
              o[2] := -0.0045942820899910*o[1];
              o[3] := 0.00217746787145710 + o[2];
              o[4] := o[3]*g.tau;
              o[5] := o[1]*g.tau;
              o[6] := o[1]*o[1];
              o[7] := o[6]*o[6];
              o[8] := o[7]*g.tau;
              o[9] := -7.9449656719138e-6*o[8];
              o[10] := g.pi*g.pi;
              o[11] := -0.0137828462699730*o[1];
              g.g := g.pi*(-0.000125631835895920 + o[4] + g.pi*(-3.9724828359569e-6*o[8] + 1.29192282897840e-7*o[5]*g.pi)) + (-0.0248051489334660 + g.tau*(0.36901534980333 + g.tau*(-3.11613182139250 + g.tau*(-13.1799836742010 + (6.8540841634434 - 0.32961626538917*g.tau)*g.tau + Modelica.Math.log(g.pi)))))/o[5];
              g.gpi := (1.0 + g.pi*(-0.000125631835895920 + o[4] + g.pi*(o[9] + 3.8757684869352e-7*o[5]*g.pi)))/g.pi;
              g.gpipi := (-1.00000000000000 + o[10]*(o[9] + 7.7515369738704e-7*o[5]*g.pi))/o[10];
              g.gtau := g.pi*(0.00217746787145710 + o[11] + g.pi*(-0.000035752345523612*o[7] + 3.8757684869352e-7*o[1]*g.pi)) + (0.074415446800398 + g.tau*(-0.73803069960666 + (3.11613182139250 + o[1]*(6.8540841634434 - 0.65923253077834*g.tau))*g.tau))/o[6];
              g.gtautau := (-0.297661787201592 + g.tau*(2.21409209881998 + (-6.2322636427850 - 0.65923253077834*o[5])*g.tau))/(o[6]*g.tau) + g.pi*(-0.0275656925399460*g.tau + g.pi*(-0.000286018764188897*o[1]*o[6]*g.tau + 7.7515369738704e-7*g.pi*g.tau));
              g.gtaupi := 0.00217746787145710 + o[11] + g.pi*(-0.000071504691047224*o[7] + 1.16273054608056e-6*o[1]*g.pi);
            end g5;

            function tph1
              input SI.Pressure p;
              input SI.SpecificEnthalpy h;
              output SI.Temperature T;
            protected
              Real pi;
              Real eta1;
              Real[3] o;
            algorithm
              assert(p > triple.ptriple, "IF97 medium function tph1 called with too low pressure\n" + "p = " + String(p) + " Pa <= " + String(triple.ptriple) + " Pa (triple point pressure)");
              pi := p/data.PSTAR2;
              eta1 := h/data.HSTAR1 + 1.0;
              o[1] := eta1*eta1;
              o[2] := o[1]*o[1];
              o[3] := o[2]*o[2];
              T := -238.724899245210 - 13.3917448726020*pi + eta1*(404.21188637945 + 43.211039183559*pi + eta1*(113.497468817180 - 54.010067170506*pi + eta1*(30.5358922039160*pi + eta1*(-6.5964749423638*pi + o[1]*(-5.8457616048039 + o[2]*(pi*(0.0093965400878363 + (-0.0000258586412820730 + 6.6456186191635e-8*pi)*pi) + o[2]*o[3]*(-0.000152854824131400 + o[1]*o[3]*(-1.08667076953770e-6 + pi*(1.15736475053400e-7 + pi*(-4.0644363084799e-9 + pi*(8.0670734103027e-11 + pi*(-9.3477771213947e-13 + (5.8265442020601e-15 - 1.50201859535030e-17*pi)*pi))))))))))));
            end tph1;

            function tph2
              input SI.Pressure p;
              input SI.SpecificEnthalpy h;
              output SI.Temperature T;
            protected
              Real pi;
              Real pi2b;
              Real pi2c;
              Real eta;
              Real etabc;
              Real eta2a;
              Real eta2b;
              Real eta2c;
              Real[8] o;
            algorithm
              pi := p*data.IPSTAR;
              eta := h*data.IHSTAR;
              etabc := h*1.0e-3;
              if (pi < 4.0) then
                eta2a := eta - 2.1;
                o[1] := eta2a*eta2a;
                o[2] := o[1]*o[1];
                o[3] := pi*pi;
                o[4] := o[3]*o[3];
                o[5] := o[3]*pi;
                T := 1089.89523182880 + (1.84457493557900 - 0.0061707422868339*pi)*pi + eta2a*(849.51654495535 - 4.1792700549624*pi + eta2a*(-107.817480918260 + (6.2478196935812 - 0.310780466295830*pi)*pi + eta2a*(33.153654801263 - 17.3445631081140*pi + o[2]*(-7.4232016790248 + pi*(-200.581768620960 + 11.6708730771070*pi) + o[1]*(271.960654737960*pi + o[1]*(-455.11318285818*pi + eta2a*(1.38657242832260*o[4] + o[1]*o[2]*(3091.96886047550*pi + o[1]*(11.7650487243560 + o[2]*(-13551.3342407750*o[5] + o[2]*(-62.459855192507*o[3]*o[4]*pi + o[2]*(o[4]*(235988.325565140 + 7399.9835474766*pi) + o[1]*(19127.7292396600*o[3]*o[4] + o[1]*(o[3]*(1.28127984040460e8 - 551966.97030060*o[5]) + o[1]*(-9.8554909623276e8*o[3] + o[1]*(2.82245469730020e9*o[3] + o[1]*(o[3]*(-3.5948971410703e9 + 3.7154085996233e6*o[5]) + o[1]*pi*(252266.403578720 + pi*(1.72273499131970e9 + pi*(1.28487346646500e7 + (-1.31052365450540e7 - 415351.64835634*o[3])*pi))))))))))))))))))));
              elseif (pi < (0.12809002730136e-03*etabc - 0.67955786399241)*etabc + 0.90584278514723e3) then
                eta2b := eta - 2.6;
                pi2b := pi - 2.0;
                o[1] := pi2b*pi2b;
                o[2] := o[1]*pi2b;
                o[3] := o[1]*o[1];
                o[4] := eta2b*eta2b;
                o[5] := o[4]*o[4];
                o[6] := o[4]*o[5];
                o[7] := o[5]*o[5];
                T := 1489.50410795160 + 0.93747147377932*pi2b + eta2b*(743.07798314034 + o[2]*(0.000110328317899990 - 1.75652339694070e-18*o[1]*o[3]) + eta2b*(-97.708318797837 + pi2b*(3.3593118604916 + pi2b*(-0.0218107553247610 + pi2b*(0.000189552483879020 + (2.86402374774560e-7 - 8.1456365207833e-14*o[2])*pi2b))) + o[5]*(3.3809355601454*pi2b + o[4]*(-0.108297844036770*o[1] + o[5]*(2.47424647056740 + (0.168445396719040 + o[1]*(0.00308915411605370 - 0.0000107798573575120*pi2b))*pi2b + o[6]*(-0.63281320016026 + pi2b*(0.73875745236695 + (-0.046333324635812 + o[1]*(-0.000076462712454814 + 2.82172816350400e-7*pi2b))*pi2b) + o[6]*(1.13859521296580 + pi2b*(-0.47128737436186 + o[1]*(0.00135555045549490 + (0.0000140523928183160 + 1.27049022719450e-6*pi2b)*pi2b)) + o[5]*(-0.47811863648625 + (0.150202731397070 + o[2]*(-0.0000310838143314340 + o[1]*(-1.10301392389090e-8 - 2.51805456829620e-11*pi2b)))*pi2b + o[5]*o[7]*(0.0085208123431544 + pi2b*(-0.00217641142197500 + pi2b*(0.000071280351959551 + o[1]*(-1.03027382121030e-6 + (7.3803353468292e-8 + 8.6934156344163e-15*o[3])*pi2b))))))))))));
              else
                eta2c := eta - 1.8;
                pi2c := pi + 25.0;
                o[1] := pi2c*pi2c;
                o[2] := o[1]*o[1];
                o[3] := o[1]*o[2]*pi2c;
                o[4] := 1/o[3];
                o[5] := o[1]*o[2];
                o[6] := eta2c*eta2c;
                o[7] := o[2]*o[2];
                o[8] := o[6]*o[6];
                T := eta2c*((859777.22535580 + o[1]*(482.19755109255 + 1.12615974072300e-12*o[5]))/o[1] + eta2c*((-5.8340131851590e11 + (2.08255445631710e10 + 31081.0884227140*o[2])*pi2c)/o[5] + o[6]*(o[8]*(o[6]*(1.23245796908320e-7*o[5] + o[6]*(-1.16069211309840e-6*o[5] + o[8]*(0.0000278463670885540*o[5] + (-0.00059270038474176*o[5] + 0.00129185829918780*o[5]*o[6])*o[8]))) - 10.8429848800770*pi2c) + o[4]*(7.3263350902181e12 + o[7]*(3.7966001272486 + (-0.045364172676660 - 1.78049822406860e-11*o[2])*pi2c))))) + o[4]*(-3.2368398555242e12 + pi2c*(3.5825089945447e11 + pi2c*(-1.07830682174700e10 + o[1]*pi2c*(610747.83564516 + pi2c*(-25745.7236041700 + (1208.23158659360 + 1.45591156586980e-13*o[5])*pi2c)))));
              end if;
            end tph2;

            function tsat
              input SI.Pressure p;
              output SI.Temperature t_sat;
            protected
              Real pi;
              Real[20] o;
            algorithm
              assert(p > triple.ptriple, "IF97 medium function tsat called with too low pressure\n" + "p = " + String(p) + " Pa <= " + String(triple.ptriple) + " Pa (triple point pressure)");
              pi := min(p, data.PCRIT)*data.IPSTAR;
              o[1] := pi^0.25;
              o[2] := -3.2325550322333e6*o[1];
              o[3] := pi^0.5;
              o[4] := -724213.16703206*o[3];
              o[5] := 405113.40542057 + o[2] + o[4];
              o[6] := -17.0738469400920*o[1];
              o[7] := 14.9151086135300 + o[3] + o[6];
              o[8] := -4.0*o[5]*o[7];
              o[9] := 12020.8247024700*o[1];
              o[10] := 1167.05214527670*o[3];
              o[11] := -4823.2657361591 + o[10] + o[9];
              o[12] := o[11]*o[11];
              o[13] := o[12] + o[8];
              o[14] := o[13]^0.5;
              o[15] := -o[14];
              o[16] := -12020.8247024700*o[1];
              o[17] := -1167.05214527670*o[3];
              o[18] := 4823.2657361591 + o[15] + o[16] + o[17];
              o[19] := 1/o[18];
              o[20] := 2.0*o[19]*o[5];
              t_sat := 0.5*(650.17534844798 + o[20] - (-4.0*(-0.238555575678490 + 1300.35069689596*o[19]*o[5]) + (650.17534844798 + o[20])^2.0)^0.5);
              annotation(derivative = tsat_der);
            end tsat;

            function dtsatofp
              input SI.Pressure p;
              output Real dtsat(unit = "K/Pa");
            protected
              Real pi;
              Real[49] o;
            algorithm
              pi := max(Modelica.Constants.small, p*data.IPSTAR);
              o[1] := pi^0.75;
              o[2] := 1/o[1];
              o[3] := -4.268461735023*o[2];
              o[4] := sqrt(pi);
              o[5] := 1/o[4];
              o[6] := 0.5*o[5];
              o[7] := o[3] + o[6];
              o[8] := pi^0.25;
              o[9] := -3.2325550322333e6*o[8];
              o[10] := -724213.16703206*o[4];
              o[11] := 405113.40542057 + o[10] + o[9];
              o[12] := -4*o[11]*o[7];
              o[13] := -808138.758058325*o[2];
              o[14] := -362106.58351603*o[5];
              o[15] := o[13] + o[14];
              o[16] := -17.073846940092*o[8];
              o[17] := 14.91510861353 + o[16] + o[4];
              o[18] := -4*o[15]*o[17];
              o[19] := 3005.2061756175*o[2];
              o[20] := 583.52607263835*o[5];
              o[21] := o[19] + o[20];
              o[22] := 12020.82470247*o[8];
              o[23] := 1167.0521452767*o[4];
              o[24] := -4823.2657361591 + o[22] + o[23];
              o[25] := 2.0*o[21]*o[24];
              o[26] := o[12] + o[18] + o[25];
              o[27] := -4.0*o[11]*o[17];
              o[28] := o[24]*o[24];
              o[29] := o[27] + o[28];
              o[30] := sqrt(o[29]);
              o[31] := 1/o[30];
              o[32] := (-o[30]);
              o[33] := -12020.82470247*o[8];
              o[34] := -1167.0521452767*o[4];
              o[35] := 4823.2657361591 + o[32] + o[33] + o[34];
              o[36] := o[30];
              o[37] := -4823.2657361591 + o[22] + o[23] + o[36];
              o[38] := o[37]*o[37];
              o[39] := 1/o[38];
              o[40] := -1.72207339365771*o[30];
              o[41] := 21592.2055343628*o[8];
              o[42] := o[30]*o[8];
              o[43] := -8192.87114842946*o[4];
              o[44] := -0.510632954559659*o[30]*o[4];
              o[45] := -3100.02526152368*o[1];
              o[46] := pi;
              o[47] := 1295.95640782102*o[46];
              o[48] := 2862.09212505088 + o[40] + o[41] + o[42] + o[43] + o[44] + o[45] + o[47];
              o[49] := 1/(o[35]*o[35]);
              dtsat := data.IPSTAR*0.5*((2.0*o[15])/o[35] - 2.*o[11]*(-3005.2061756175*o[2] - 0.5*o[26]*o[31] - 583.52607263835*o[5])*o[49] - (20953.46356643991*(o[39]*(1295.95640782102 + 5398.05138359071*o[2] + 0.25*o[2]*o[30] - 0.861036696828853*o[26]*o[31] - 0.255316477279829*o[26]*o[31]*o[4] - 4096.43557421473*o[5] - 0.255316477279829*o[30]*o[5] - 2325.01894614276/o[8] + 0.5*o[26]*o[31]*o[8]) - 2.0*(o[19] + o[20] + 0.5*o[26]*o[31])*o[48]*o[37]^(-3)))/sqrt(o[39]*o[48]));
            end dtsatofp;

            function tsat_der
              input SI.Pressure p;
              input Real der_p(unit = "Pa/s");
              output Real der_tsat(unit = "K/s");
            protected
              Real dtp;
            algorithm
              dtp := dtsatofp(p);
              der_tsat := dtp*der_p;
            end tsat_der;

            function psat
              input SI.Temperature T;
              output SI.Pressure p_sat;
            protected
              Real[8] o;
              Real Tlim = min(T, data.TCRIT);
            algorithm
              assert(T >= 273.16, "IF97 medium function psat: input temperature (= " + String(triple.Ttriple) + " K).\n" + "lower than the triple point temperature 273.16 K");
              o[1] := -650.17534844798 + Tlim;
              o[2] := 1/o[1];
              o[3] := -0.238555575678490*o[2];
              o[4] := o[3] + Tlim;
              o[5] := -4823.2657361591*o[4];
              o[6] := o[4]*o[4];
              o[7] := 14.9151086135300*o[6];
              o[8] := 405113.40542057 + o[5] + o[7];
              p_sat := 16.0e6*o[8]*o[8]*o[8]*o[8]*1/(3.2325550322333e6 - 12020.8247024700*o[4] + 17.0738469400920*o[6] + (-4.0*(-724213.16703206 + 1167.05214527670*o[4] + o[6])*o[8] + (-3.2325550322333e6 + 12020.8247024700*o[4] - 17.0738469400920*o[6])^2.0)^0.5)^4.0;
              annotation(derivative = psat_der);
            end psat;

            function dptofT
              input SI.Temperature T;
              output Real dpt(unit = "Pa/K");
            protected
              Real[31] o;
              Real Tlim;
            algorithm
              Tlim := min(T, data.TCRIT);
              o[1] := -650.17534844798 + Tlim;
              o[2] := 1/o[1];
              o[3] := -0.238555575678490*o[2];
              o[4] := o[3] + Tlim;
              o[5] := -4823.2657361591*o[4];
              o[6] := o[4]*o[4];
              o[7] := 14.9151086135300*o[6];
              o[8] := 405113.40542057 + o[5] + o[7];
              o[9] := o[8]*o[8];
              o[10] := o[9]*o[9];
              o[11] := o[1]*o[1];
              o[12] := 1/o[11];
              o[13] := 0.238555575678490*o[12];
              o[14] := 1.00000000000000 + o[13];
              o[15] := 12020.8247024700*o[4];
              o[16] := -17.0738469400920*o[6];
              o[17] := -3.2325550322333e6 + o[15] + o[16];
              o[18] := -4823.2657361591*o[14];
              o[19] := 29.8302172270600*o[14]*o[4];
              o[20] := o[18] + o[19];
              o[21] := 1167.05214527670*o[4];
              o[22] := -724213.16703206 + o[21] + o[6];
              o[23] := o[17]*o[17];
              o[24] := -4.0000000000000*o[22]*o[8];
              o[25] := o[23] + o[24];
              o[26] := sqrt(o[25]);
              o[27] := -12020.8247024700*o[4];
              o[28] := 17.0738469400920*o[6];
              o[29] := 3.2325550322333e6 + o[26] + o[27] + o[28];
              o[30] := o[29]*o[29];
              o[31] := o[30]*o[30];
              dpt := 1e6*((-64.0*o[10]*(-12020.8247024700*o[14] + 34.147693880184*o[14]*o[4] + (0.5*(-4.0*o[20]*o[22] + 2.00000000000000*o[17]*(12020.8247024700*o[14] - 34.147693880184*o[14]*o[4]) - 4.0*(1167.05214527670*o[14] + 2.0*o[14]*o[4])*o[8]))/o[26]))/(o[29]*o[31]) + (64.*o[20]*o[8]*o[9])/o[31]);
            end dptofT;

            function psat_der
              input SI.Temperature T;
              input Real der_T(unit = "K/s");
              output Real der_psat(unit = "Pa/s");
            protected
              Real dpt;
            algorithm
              dpt := dptofT(T);
              der_psat := dpt*der_T;
            end psat_der;

            function h3ab_p
              output SI.SpecificEnthalpy h;
              input SI.Pressure p;
            protected
              constant Real[:] n = {0.201464004206875e4, 0.374696550136983e1, -0.219921901054187e-1, 0.875131686009950e-4};
              constant SI.SpecificEnthalpy hstar = 1000;
              constant SI.Pressure pstar = 1e6;
              Real pi = p/pstar;
            algorithm
              h := (n[1] + n[2]*pi + n[3]*pi^2 + n[4]*pi^3)*hstar;
            end h3ab_p;

            function T3a_ph
              input SI.Pressure p;
              input SI.SpecificEnthalpy h;
              output SI.Temperature T;
            protected
              constant Real[:] n = {-0.133645667811215e-6, 0.455912656802978e-5, -0.146294640700979e-4, 0.639341312970080e-2, 0.372783927268847e3, -0.718654377460447e4, 0.573494752103400e6, -0.267569329111439e7, -0.334066283302614e-4, -0.245479214069597e-1, 0.478087847764996e2, 0.764664131818904e-5, 0.128350627676972e-2, 0.171219081377331e-1, -0.851007304583213e1, -0.136513461629781e-1, -0.384460997596657e-5, 0.337423807911655e-2, -0.551624873066791, 0.729202277107470, -0.992522757376041e-2, -0.119308831407288, 0.793929190615421, 0.454270731799386, 0.209998591259910, -0.642109823904738e-2, -0.235155868604540e-1, 0.252233108341612e-2, -0.764885133368119e-2, 0.136176427574291e-1, -0.133027883575669e-1};
              constant Real[:] I = {-12, -12, -12, -12, -12, -12, -12, -12, -10, -10, -10, -8, -8, -8, -8, -5, -3, -2, -2, -2, -1, -1, 0, 0, 1, 3, 3, 4, 4, 10, 12};
              constant Real[:] J = {0, 1, 2, 6, 14, 16, 20, 22, 1, 5, 12, 0, 2, 4, 10, 2, 0, 1, 3, 4, 0, 2, 0, 1, 1, 0, 1, 0, 3, 4, 5};
              constant SI.SpecificEnthalpy hstar = 2300e3;
              constant SI.Pressure pstar = 100e6;
              constant SI.Temperature Tstar = 760;
              Real pi = p/pstar;
              Real eta = h/hstar;
            algorithm
              T := sum(n[i]*(pi + 0.240)^I[i]*(eta - 0.615)^J[i] for i in 1:31)*Tstar;
            end T3a_ph;

            function T3b_ph
              input SI.Pressure p;
              input SI.SpecificEnthalpy h;
              output SI.Temperature T;
            protected
              constant Real[:] n = {0.323254573644920e-4, -0.127575556587181e-3, -0.475851877356068e-3, 0.156183014181602e-2, 0.105724860113781, -0.858514221132534e2, 0.724140095480911e3, 0.296475810273257e-2, -0.592721983365988e-2, -0.126305422818666e-1, -0.115716196364853, 0.849000969739595e2, -0.108602260086615e-1, 0.154304475328851e-1, 0.750455441524466e-1, 0.252520973612982e-1, -0.602507901232996e-1, -0.307622221350501e1, -0.574011959864879e-1, 0.503471360939849e1, -0.925081888584834, 0.391733882917546e1, -0.773146007130190e2, 0.949308762098587e4, -0.141043719679409e7, 0.849166230819026e7, 0.861095729446704, 0.323346442811720, 0.873281936020439, -0.436653048526683, 0.286596714529479, -0.131778331276228, 0.676682064330275e-2};
              constant Real[:] I = {-12, -12, -10, -10, -10, -10, -10, -8, -8, -8, -8, -8, -6, -6, -6, -4, -4, -3, -2, -2, -1, -1, -1, -1, -1, -1, 0, 0, 1, 3, 5, 6, 8};
              constant Real[:] J = {0, 1, 0, 1, 5, 10, 12, 0, 1, 2, 4, 10, 0, 1, 2, 0, 1, 5, 0, 4, 2, 4, 6, 10, 14, 16, 0, 2, 1, 1, 1, 1, 1};
              constant SI.Temperature Tstar = 860;
              constant SI.Pressure pstar = 100e6;
              constant SI.SpecificEnthalpy hstar = 2800e3;
              Real pi = p/pstar;
              Real eta = h/hstar;
            algorithm
              T := sum(n[i]*(pi + 0.298)^I[i]*(eta - 0.720)^J[i] for i in 1:33)*Tstar;
            end T3b_ph;

            function v3a_ph
              input SI.Pressure p;
              input SI.SpecificEnthalpy h;
              output SI.SpecificVolume v;
            protected
              constant Real[:] n = {0.529944062966028e-2, -0.170099690234461, 0.111323814312927e2, -0.217898123145125e4, -0.506061827980875e-3, 0.556495239685324, -0.943672726094016e1, -0.297856807561527, 0.939353943717186e2, 0.192944939465981e-1, 0.421740664704763, -0.368914126282330e7, -0.737566847600639e-2, -0.354753242424366, -0.199768169338727e1, 0.115456297059049e1, 0.568366875815960e4, 0.808169540124668e-2, 0.172416341519307, 0.104270175292927e1, -0.297691372792847, 0.560394465163593, 0.275234661176914, -0.148347894866012, -0.651142513478515e-1, -0.292468715386302e1, 0.664876096952665e-1, 0.352335014263844e1, -0.146340792313332e-1, -0.224503486668184e1, 0.110533464706142e1, -0.408757344495612e-1};
              constant Real[:] I = {-12, -12, -12, -12, -10, -10, -10, -8, -8, -6, -6, -6, -4, -4, -3, -2, -2, -1, -1, -1, -1, 0, 0, 1, 1, 1, 2, 2, 3, 4, 5, 8};
              constant Real[:] J = {6, 8, 12, 18, 4, 7, 10, 5, 12, 3, 4, 22, 2, 3, 7, 3, 16, 0, 1, 2, 3, 0, 1, 0, 1, 2, 0, 2, 0, 2, 2, 2};
              constant SI.Volume vstar = 0.0028;
              constant SI.Pressure pstar = 100e6;
              constant SI.SpecificEnthalpy hstar = 2100e3;
              Real pi = p/pstar;
              Real eta = h/hstar;
            algorithm
              v := sum(n[i]*(pi + 0.128)^I[i]*(eta - 0.727)^J[i] for i in 1:32)*vstar;
            end v3a_ph;

            function v3b_ph
              input SI.Pressure p;
              input SI.SpecificEnthalpy h;
              output SI.SpecificVolume v;
            protected
              constant Real[:] n = {-0.225196934336318e-8, 0.140674363313486e-7, 0.233784085280560e-5, -0.331833715229001e-4, 0.107956778514318e-2, -0.271382067378863, 0.107202262490333e1, -0.853821329075382, -0.215214194340526e-4, 0.769656088222730e-3, -0.431136580433864e-2, 0.453342167309331, -0.507749535873652, -0.100475154528389e3, -0.219201924648793, -0.321087965668917e1, 0.607567815637771e3, 0.557686450685932e-3, 0.187499040029550, 0.905368030448107e-2, 0.285417173048685, 0.329924030996098e-1, 0.239897419685483, 0.482754995951394e1, -0.118035753702231e2, 0.169490044091791, -0.179967222507787e-1, 0.371810116332674e-1, -0.536288335065096e-1, 0.160697101092520e1};
              constant Real[:] I = {-12, -12, -8, -8, -8, -8, -8, -8, -6, -6, -6, -6, -6, -6, -4, -4, -4, -3, -3, -2, -2, -1, -1, -1, -1, 0, 1, 1, 2, 2};
              constant Real[:] J = {0, 1, 0, 1, 3, 6, 7, 8, 0, 1, 2, 5, 6, 10, 3, 6, 10, 0, 2, 1, 2, 0, 1, 4, 5, 0, 0, 1, 2, 6};
              constant SI.Volume vstar = 0.0088;
              constant SI.Pressure pstar = 100e6;
              constant SI.SpecificEnthalpy hstar = 2800e3;
              Real pi = p/pstar;
              Real eta = h/hstar;
            algorithm
              v := sum(n[i]*(pi + 0.0661)^I[i]*(eta - 0.720)^J[i] for i in 1:30)*vstar;
            end v3b_ph;
          end Basic;

          package Isentropic
            function hofpT1
              input SI.Pressure p;
              input SI.Temperature T;
              output SI.SpecificEnthalpy h;
            protected
              Real[13] o;
              Real pi1;
              Real tau;
              Real tau1;
            algorithm
              tau := data.TSTAR1/T;
              pi1 := 7.1 - p/data.PSTAR1;
              assert(p > triple.ptriple, "IF97 medium function hofpT1 called with too low pressure\n" + "p = " + String(p) + " Pa <= " + String(triple.ptriple) + " Pa (triple point pressure)");
              tau1 := -1.222 + tau;
              o[1] := tau1*tau1;
              o[2] := o[1]*tau1;
              o[3] := o[1]*o[1];
              o[4] := o[3]*o[3];
              o[5] := o[1]*o[4];
              o[6] := o[1]*o[3];
              o[7] := o[3]*tau1;
              o[8] := o[3]*o[4];
              o[9] := pi1*pi1;
              o[10] := o[9]*o[9];
              o[11] := o[10]*o[10];
              o[12] := o[4]*o[4];
              o[13] := o[12]*o[12];
              h := data.RH2O*T*tau*(pi1*((-0.00254871721114236 + o[1]*(0.00424944110961118 + (0.018990068218419 + (-0.021841717175414 - 0.00015851507390979*o[1])*o[1])*o[6]))/o[5] + pi1*((0.00141552963219801 + o[3]*(0.000047661393906987 + o[1]*(-0.0000132425535992538 - 1.2358149370591e-14*o[1]*o[3]*o[4])))/o[3] + pi1*((0.000126718579380216 - 5.11230768720618e-9*o[5])/o[7] + pi1*((0.000011212640954 + o[2]*(1.30342445791202e-6 - 1.4341729937924e-12*o[8]))/o[6] + pi1*(o[9]*pi1*((1.40077319158051e-8 + 1.04549227383804e-9*o[7])/o[8] + o[10]*o[11]*pi1*(1.9941018075704e-17/(o[1]*o[12]*o[3]*o[4]) + o[9]*(-4.48827542684151e-19/o[13] + o[10]*o[9]*(pi1*(4.65957282962769e-22/(o[13]*o[4]) + pi1*((3.83502057899078e-24*pi1)/(o[1]*o[13]*o[4]) - 7.2912378325616e-23/(o[13]*o[4]*tau1))) - 1.00075970318621e-21/(o[1]*o[13]*o[3]*tau1))))) + 3.24135974880936e-6/(o[4]*tau1)))))) + (-0.29265942426334 + tau1*(0.84548187169114 + o[1]*(3.3855169168385 + tau1*(-1.91583926775744 + tau1*(0.47316115539684 + (-0.066465668798004 + 0.0040607314991784*tau1)*tau1)))))/o[2]);
            end hofpT1;

            function hofpT2
              input SI.Pressure p;
              input SI.Temperature T;
              output SI.SpecificEnthalpy h;
            protected
              Real[16] o;
              Real pi;
              Real tau;
              Real tau2;
            algorithm
              assert(p > triple.ptriple, "IF97 medium function hofpT2 called with too low pressure\n" + "p = " + String(p) + " Pa <= " + String(triple.ptriple) + " Pa (triple point pressure)");
              pi := p/data.PSTAR2;
              tau := data.TSTAR2/T;
              tau2 := -0.5 + tau;
              o[1] := tau*tau;
              o[2] := o[1]*o[1];
              o[3] := tau2*tau2;
              o[4] := o[3]*tau2;
              o[5] := o[3]*o[3];
              o[6] := o[5]*o[5];
              o[7] := o[6]*o[6];
              o[8] := o[5]*o[6]*o[7]*tau2;
              o[9] := o[3]*o[5];
              o[10] := o[5]*o[6]*tau2;
              o[11] := o[3]*o[7]*tau2;
              o[12] := o[3]*o[5]*o[6];
              o[13] := o[5]*o[6]*o[7];
              o[14] := pi*pi;
              o[15] := o[14]*o[14];
              o[16] := o[7]*o[7];
              h := data.RH2O*T*tau*((0.0280439559151 + tau*(-0.2858109552582 + tau*(1.2213149471784 + tau*(-2.848163942888 + tau*(4.38395111945 + o[1]*(10.08665568018 + (-0.5681726521544 + 0.06380539059921*tau)*tau))))))/(o[1]*o[2]) + pi*(-0.017834862292358 + tau2*(-0.09199202739273 + (-0.172743777250296 - 0.30195167236758*o[4])*tau2) + pi*(-0.000033032641670203 + (-0.0003789797503263 + o[3]*(-0.015757110897342 + o[4]*(-0.306581069554011 - 0.000960283724907132*o[8])))*tau2 + pi*(4.3870667284435e-7 + o[3]*(-0.00009683303171571 + o[4]*(-0.0090203547252888 - 1.42338887469272*o[8])) + pi*(-7.8847309559367e-10 + (2.558143570457e-8 + 1.44676118155521e-6*tau2)*tau2 + pi*(0.0000160454534363627*o[9] + pi*((-5.0144299353183e-11 + o[10]*(-0.033874355714168 - 836.35096769364*o[11]))*o[3] + pi*((-0.0000138839897890111 - 0.973671060893475*o[12])*o[3]*o[6] + pi*((9.0049690883672e-11 - 296.320827232793*o[13])*o[3]*o[5]*tau2 + pi*(2.57526266427144e-7*o[5]*o[6] + pi*(o[4]*(4.1627860840696e-19 + (-1.0234747095929e-12 - 1.40254511313154e-8*o[5])*o[9]) + o[14]*o[15]*(o[13]*(-2.34560435076256e-9 + 5.3465159397045*o[5]*o[7]*tau2) + o[14]*(-19.1874828272775*o[16]*o[6]*o[7] + o[14]*(o[11]*(1.78371690710842e-23 + (1.07202609066812e-11 - 0.000201611844951398*o[10])*o[3]*o[5]*o[6]*tau2) + pi*(-1.24017662339842e-24*o[5]*o[7] + pi*(0.000200482822351322*o[16]*o[5]*o[7] + pi*(-4.97975748452559e-14*o[16]*o[3]*o[5] + o[6]*o[7]*(1.90027787547159e-27 + o[12]*(2.21658861403112e-15 - 0.0000547344301999018*o[3]*o[7]))*pi*tau2)))))))))))))))));
            end hofpT2;
          end Isentropic;

          package Inverses
            function fixdT
              input SI.Density din;
              input SI.Temperature Tin;
              output SI.Density dout;
              output SI.Temperature Tout;
            protected
              SI.Temperature Tmin;
              SI.Temperature Tmax;
            algorithm
              if (din > 765.0) then
                dout := 765.0;
              elseif (din < 110.0) then
                dout := 110.0;
              else
                dout := din;
              end if;
              if (dout < 390.0) then
                Tmax := 554.3557377 + dout*0.809344262;
              else
                Tmax := 1116.85 - dout*0.632948717;
              end if;
              if (dout < data.DCRIT) then
                Tmin := data.TCRIT*(1.0 - (dout - data.DCRIT)*(dout - data.DCRIT)/1.0e6);
              else
                Tmin := data.TCRIT*(1.0 - (dout - data.DCRIT)*(dout - data.DCRIT)/1.44e6);
              end if;
              if (Tin < Tmin) then
                Tout := Tmin;
              elseif (Tin > Tmax) then
                Tout := Tmax;
              else
                Tout := Tin;
              end if;
            end fixdT;

            function dofp13
              input SI.Pressure p;
              output SI.Density d;
            protected
              Real p2;
              Real[3] o;
            algorithm
              p2 := 7.1 - 6.04960677555959e-8*p;
              o[1] := p2*p2;
              o[2] := o[1]*o[1];
              o[3] := o[2]*o[2];
              d := 57.4756752485113/(0.0737412153522555 + p2*(0.00145092247736023 + p2*(0.000102697173772229 + p2*(0.0000114683182476084 + p2*(1.99080616601101e-6 + o[1]*p2*(1.13217858826367e-8 + o[2]*o[3]*p2*(1.35549330686006e-17 + o[1]*(-3.11228834832975e-19 + o[1]*o[2]*(-7.02987180039442e-22 + p2*(3.29199117056433e-22 + (-5.17859076694812e-23 + 2.73712834080283e-24*p2)*p2))))))))));
            end dofp13;

            function dofp23
              input SI.Pressure p;
              output SI.Density d;
            protected
              SI.Temperature T;
              Real[13] o;
              Real taug;
              Real pi;
              Real gpi23;
            algorithm
              pi := p/data.PSTAR2;
              T := 572.54459862746 + 31.3220101646784*(-13.91883977887 + pi)^0.5;
              o[1] := (-13.91883977887 + pi)^0.5;
              taug := -0.5 + 540.0/(572.54459862746 + 31.3220101646784*o[1]);
              o[2] := taug*taug;
              o[3] := o[2]*taug;
              o[4] := o[2]*o[2];
              o[5] := o[4]*o[4];
              o[6] := o[5]*o[5];
              o[7] := o[4]*o[5]*o[6]*taug;
              o[8] := o[4]*o[5]*taug;
              o[9] := o[2]*o[4]*o[5];
              o[10] := pi*pi;
              o[11] := o[10]*o[10];
              o[12] := o[4]*o[6]*taug;
              o[13] := o[6]*o[6];
              gpi23 := (1.0 + pi*(-0.0017731742473213 + taug*(-0.017834862292358 + taug*(-0.045996013696365 + (-0.057581259083432 - 0.05032527872793*o[3])*taug)) + pi*(taug*(-0.000066065283340406 + (-0.0003789797503263 + o[2]*(-0.007878555448671 + o[3]*(-0.087594591301146 - 0.000053349095828174*o[7])))*taug) + pi*(6.1445213076927e-8 + (1.31612001853305e-6 + o[2]*(-0.00009683303171571 + o[3]*(-0.0045101773626444 - 0.122004760687947*o[7])))*taug + pi*(taug*(-3.15389238237468e-9 + (5.116287140914e-8 + 1.92901490874028e-6*taug)*taug) + pi*(0.0000114610381688305*o[2]*o[4]*taug + pi*(o[3]*(-1.00288598706366e-10 + o[8]*(-0.012702883392813 - 143.374451604624*o[2]*o[6]*taug)) + pi*(-4.1341695026989e-17 + o[2]*o[5]*(-8.8352662293707e-6 - 0.272627897050173*o[9])*taug + pi*(o[5]*(9.0049690883672e-11 - 65.8490727183984*o[4]*o[5]*o[6]) + pi*(1.78287415218792e-7*o[8] + pi*(o[4]*(1.0406965210174e-18 + o[2]*(-1.0234747095929e-12 - 1.0018179379511e-8*o[4])*o[4]) + o[10]*o[11]*((-1.29412653835176e-9 + 1.71088510070544*o[12])*o[7] + o[10]*(-6.05920510335078*o[13]*o[5]*o[6]*taug + o[10]*(o[4]*o[6]*(1.78371690710842e-23 + o[2]*o[4]*o[5]*(6.1258633752464e-12 - 0.000084004935396416*o[8])*taug) + pi*(-1.24017662339842e-24*o[12] + pi*(0.0000832192847496054*o[13]*o[4]*o[6]*taug + pi*(o[2]*o[5]*o[6]*(1.75410265428146e-27 + (1.32995316841867e-15 - 0.0000226487297378904*o[2]*o[6])*o[9])*pi - 2.93678005497663e-14*o[13]*o[2]*o[4]*taug)))))))))))))))))/pi;
              d := p/(data.RH2O*T*pi*gpi23);
            end dofp23;

            function dofpt3
              input SI.Pressure p;
              input SI.Temperature T;
              input SI.Pressure delp;
              output SI.Density d;
              output Integer error = 0;
            protected
              SI.Density dguess;
              Integer i = 0;
              Real dp;
              SI.Density deld;
              Modelica.Media.Common.HelmholtzDerivs f;
              Modelica.Media.Common.NewtonDerivatives_pT nDerivs;
              Boolean found = false;
              Boolean supercritical;
              Boolean liquid;
              SI.Density dmin;
              SI.Density dmax;
              SI.Temperature Tmax;
              Real damping;
            algorithm
              found := false;
              assert(p >= data.PLIMIT4A, "BaseIF97.dofpt3: function called outside of region 3! p too low\n" + "p = " + String(p) + " Pa < " + String(data.PLIMIT4A) + " Pa");
              assert(T >= data.TLIMIT1, "BaseIF97.dofpt3: function called outside of region 3! T too low\n" + "T = " + String(T) + " K < " + String(data.TLIMIT1) + " K");
              assert(p >= Regions.boundary23ofT(T), "BaseIF97.dofpt3: function called outside of region 3! T too high\n" + "p = " + String(p) + " Pa, T = " + String(T) + " K");
              supercritical := p > data.PCRIT;
              damping := if supercritical then 1.0 else 1.0;
              Tmax := Regions.boundary23ofp(p);
              if supercritical then
                dmax := dofp13(p);
                dmin := dofp23(p);
                dguess := dmax - (T - data.TLIMIT1)/(data.TLIMIT1 - Tmax)*(dmax - dmin);
              else
                liquid := T < Basic.tsat(p);
                if liquid then
                  dmax := dofp13(p);
                  dmin := Regions.rhol_p_R4b(p);
                  dguess := 1.1*Regions.rhol_T(T);
                else
                  dmax := Regions.rhov_p_R4b(p);
                  dmin := dofp23(p);
                  dguess := 0.9*Regions.rhov_T(T);
                end if;
              end if;
              while ((i < IterationData.IMAX) and not found) loop
                d := dguess;
                f := Basic.f3(d, T);
                nDerivs := Modelica.Media.Common.Helmholtz_pT(f);
                dp := nDerivs.p - p;
                if (abs(dp/p) <= delp) then
                  found := true;
                else
                end if;
                deld := dp/nDerivs.pd*damping;
                d := d - deld;
                if d > dmin and d < dmax then
                  dguess := d;
                else
                  if d > dmax then
                    dguess := dmax - sqrt(Modelica.Constants.eps);
                  else
                    dguess := dmin + sqrt(Modelica.Constants.eps);
                  end if;
                end if;
                i := i + 1;
              end while;
              if not found then
                error := 1;
              else
              end if;
              assert(error <> 1, "Error in inverse function dofpt3: iteration failed");
            end dofpt3;

            function dtofph3
              input SI.Pressure p;
              input SI.SpecificEnthalpy h;
              input SI.Pressure delp;
              input SI.SpecificEnthalpy delh;
              output SI.Density d;
              output SI.Temperature T;
              output Integer error;
            protected
              SI.Temperature Tguess;
              SI.Density dguess;
              Integer i;
              Real dh;
              Real dp;
              Real det;
              Real deld;
              Real delt;
              Modelica.Media.Common.HelmholtzDerivs f;
              Modelica.Media.Common.NewtonDerivatives_ph nDerivs;
              Boolean found = false;
              Integer subregion;
            algorithm
              if p < data.PCRIT then
                subregion := if h < (Regions.hl_p(p) + 10.0) then 1 else if h > (Regions.hv_p(p) - 10.0) then 2 else 0;
                assert(subregion <> 0, "Inverse iteration of dt from ph called in 2 phase region: this can not work");
              else
                subregion := if h < Basic.h3ab_p(p) then 1 else 2;
              end if;
              T := if subregion == 1 then Basic.T3a_ph(p, h) else Basic.T3b_ph(p, h);
              d := if subregion == 1 then 1/Basic.v3a_ph(p, h) else 1/Basic.v3b_ph(p, h);
              i := 0;
              error := 0;
              while ((i < IterationData.IMAX) and not found) loop
                f := Basic.f3(d, T);
                nDerivs := Modelica.Media.Common.Helmholtz_ph(f);
                dh := nDerivs.h - h;
                dp := nDerivs.p - p;
                if ((abs(dh/h) <= delh) and (abs(dp/p) <= delp)) then
                  found := true;
                else
                end if;
                det := nDerivs.ht*nDerivs.pd - nDerivs.pt*nDerivs.hd;
                delt := (nDerivs.pd*dh - nDerivs.hd*dp)/det;
                deld := (nDerivs.ht*dp - nDerivs.pt*dh)/det;
                T := T - delt;
                d := d - deld;
                dguess := d;
                Tguess := T;
                i := i + 1;
                (d, T) := fixdT(dguess, Tguess);
              end while;
              if not found then
                error := 1;
              else
              end if;
              assert(error <> 1, "Error in inverse function dtofph3: iteration failed");
            end dtofph3;

            function pofdt125
              input SI.Density d;
              input SI.Temperature T;
              input SI.Pressure reldd;
              input Integer region;
              output SI.Pressure p;
              output Integer error;
            protected
              Integer i;
              Modelica.Media.Common.GibbsDerivs g;
              Boolean found;
              Real dd;
              Real delp;
              Real relerr;
              SI.Pressure pguess1 = 1.0e6;
              SI.Pressure pguess2;
              constant SI.Pressure pguess5 = 0.5e6;
            algorithm
              i := 0;
              error := 0;
              pguess2 := 42800*d;
              found := false;
              if region == 1 then
                p := pguess1;
              elseif region == 2 then
                p := pguess2;
              else
                p := pguess5;
              end if;
              while ((i < IterationData.IMAX) and not found) loop
                if region == 1 then
                  g := Basic.g1(p, T);
                elseif region == 2 then
                  g := Basic.g2(p, T);
                else
                  g := Basic.g5(p, T);
                end if;
                dd := p/(data.RH2O*T*g.pi*g.gpi) - d;
                relerr := dd/d;
                if (abs(relerr) < reldd) then
                  found := true;
                else
                end if;
                delp := dd*(-p*p/(d*d*data.RH2O*T*g.pi*g.pi*g.gpipi));
                p := p - delp;
                i := i + 1;
                if not found then
                  if p < triple.ptriple then
                    p := 2.0*triple.ptriple;
                  else
                  end if;
                  if p > data.PLIMIT1 then
                    p := 0.95*data.PLIMIT1;
                  else
                  end if;
                else
                end if;
              end while;
              if not found then
                error := 1;
              else
              end if;
              assert(error <> 1, "Error in inverse function pofdt125: iteration failed");
            end pofdt125;

            function tofph5
              input SI.Pressure p;
              input SI.SpecificEnthalpy h;
              input SI.SpecificEnthalpy reldh;
              output SI.Temperature T;
              output Integer error;
            protected
              Modelica.Media.Common.GibbsDerivs g;
              SI.SpecificEnthalpy proh;
              constant SI.Temperature Tguess = 1500;
              Integer i;
              Real relerr;
              Real dh;
              Real dT;
              Boolean found;
            algorithm
              i := 0;
              error := 0;
              T := Tguess;
              found := false;
              while ((i < IterationData.IMAX) and not found) loop
                g := Basic.g5(p, T);
                proh := data.RH2O*T*g.tau*g.gtau;
                dh := proh - h;
                relerr := dh/h;
                if (abs(relerr) < reldh) then
                  found := true;
                else
                end if;
                dT := dh/(-data.RH2O*g.tau*g.tau*g.gtautau);
                T := T - dT;
                i := i + 1;
              end while;
              if not found then
                error := 1;
              else
              end if;
              assert(error <> 1, "Error in inverse function tofph5: iteration failed");
            end tofph5;
          end Inverses;

          package TwoPhase end TwoPhase;
        end BaseIF97;

        function waterBaseProp_ph
          input SI.Pressure p;
          input SI.SpecificEnthalpy h;
          input Integer phase = 0;
          input Integer region = 0;
          output Common.IF97BaseTwoPhase aux;
        protected
          Common.GibbsDerivs g;
          Common.HelmholtzDerivs f;
          Integer error;
          SI.SpecificEnthalpy h_liq;
          SI.Density d_liq;
          SI.SpecificEnthalpy h_vap;
          SI.Density d_vap;
          Common.PhaseBoundaryProperties liq;
          Common.PhaseBoundaryProperties vap;
          Common.GibbsDerivs gl;
          Common.GibbsDerivs gv;
          Modelica.Media.Common.HelmholtzDerivs fl;
          Modelica.Media.Common.HelmholtzDerivs fv;
          SI.Temperature t1;
          SI.Temperature t2;
        algorithm
          aux.region := if region == 0 then (if phase == 2 then 4 else BaseIF97.Regions.region_ph(p = p, h = h, phase = phase)) else region;
          aux.phase := if phase <> 0 then phase else if aux.region == 4 then 2 else 1;
          aux.p := max(p, 611.657);
          aux.h := max(h, 1e3);
          aux.R_s := BaseIF97.data.RH2O;
          aux.vt := 0.0;
          aux.vp := 0.0;
          if (aux.region == 1) then
            aux.T := BaseIF97.Basic.tph1(aux.p, aux.h);
            g := BaseIF97.Basic.g1(p, aux.T);
            aux.s := aux.R_s*(g.tau*g.gtau - g.g);
            aux.rho := p/(aux.R_s*aux.T*g.pi*g.gpi);
            aux.vt := aux.R_s/p*(g.pi*g.gpi - g.tau*g.pi*g.gtaupi);
            aux.pt := -g.p/g.T*(g.gpi - g.tau*g.gtaupi)/(g.gpipi*g.pi);
            aux.pd := -g.R_s*g.T*g.gpi*g.gpi/(g.gpipi);
            aux.vp := aux.R_s*aux.T/(p*p)*g.pi*g.pi*g.gpipi;
            aux.cp := -aux.R_s*g.tau*g.tau*g.gtautau;
            aux.cv := aux.R_s*(-g.tau*g.tau*g.gtautau + ((g.gpi - g.tau*g.gtaupi)*(g.gpi - g.tau*g.gtaupi)/g.gpipi));
            aux.x := 0.0;
            aux.dpT := -aux.vt/aux.vp;
          elseif (aux.region == 2) then
            aux.T := BaseIF97.Basic.tph2(aux.p, aux.h);
            g := BaseIF97.Basic.g2(p, aux.T);
            aux.s := aux.R_s*(g.tau*g.gtau - g.g);
            aux.rho := p/(aux.R_s*aux.T*g.pi*g.gpi);
            aux.vt := aux.R_s/p*(g.pi*g.gpi - g.tau*g.pi*g.gtaupi);
            aux.vp := aux.R_s*aux.T/(p*p)*g.pi*g.pi*g.gpipi;
            aux.pt := -g.p/g.T*(g.gpi - g.tau*g.gtaupi)/(g.gpipi*g.pi);
            aux.pd := -g.R_s*g.T*g.gpi*g.gpi/(g.gpipi);
            aux.cp := -aux.R_s*g.tau*g.tau*g.gtautau;
            aux.cv := aux.R_s*(-g.tau*g.tau*g.gtautau + ((g.gpi - g.tau*g.gtaupi)*(g.gpi - g.tau*g.gtaupi)/g.gpipi));
            aux.x := 1.0;
            aux.dpT := -aux.vt/aux.vp;
          elseif (aux.region == 3) then
            (aux.rho, aux.T, error) := BaseIF97.Inverses.dtofph3(p = aux.p, h = aux.h, delp = 1.0e-7, delh = 1.0e-6);
            f := BaseIF97.Basic.f3(aux.rho, aux.T);
            aux.h := aux.R_s*aux.T*(f.tau*f.ftau + f.delta*f.fdelta);
            aux.s := aux.R_s*(f.tau*f.ftau - f.f);
            aux.pd := aux.R_s*aux.T*f.delta*(2.0*f.fdelta + f.delta*f.fdeltadelta);
            aux.pt := aux.R_s*aux.rho*f.delta*(f.fdelta - f.tau*f.fdeltatau);
            aux.cv := abs(aux.R_s*(-f.tau*f.tau*f.ftautau));
            aux.cp := (aux.rho*aux.rho*aux.pd*aux.cv + aux.T*aux.pt*aux.pt)/(aux.rho*aux.rho*aux.pd);
            aux.x := 0.0;
            aux.dpT := aux.pt;
          elseif (aux.region == 4) then
            h_liq := hl_p(p);
            h_vap := hv_p(p);
            aux.x := if (h_vap <> h_liq) then (h - h_liq)/(h_vap - h_liq) else 1.0;
            if p < BaseIF97.data.PLIMIT4A then
              t1 := BaseIF97.Basic.tph1(aux.p, h_liq);
              t2 := BaseIF97.Basic.tph2(aux.p, h_vap);
              gl := BaseIF97.Basic.g1(aux.p, t1);
              gv := BaseIF97.Basic.g2(aux.p, t2);
              liq := Common.gibbsToBoundaryProps(gl);
              vap := Common.gibbsToBoundaryProps(gv);
              aux.T := t1 + aux.x*(t2 - t1);
            else
              aux.T := BaseIF97.Basic.tsat(aux.p);
              d_liq := rhol_T(aux.T);
              d_vap := rhov_T(aux.T);
              fl := BaseIF97.Basic.f3(d_liq, aux.T);
              fv := BaseIF97.Basic.f3(d_vap, aux.T);
              liq := Common.helmholtzToBoundaryProps(fl);
              vap := Common.helmholtzToBoundaryProps(fv);
            end if;
            aux.dpT := if (liq.d <> vap.d) then (vap.s - liq.s)*liq.d*vap.d/(liq.d - vap.d) else BaseIF97.Basic.dptofT(aux.T);
            aux.s := liq.s + aux.x*(vap.s - liq.s);
            aux.rho := liq.d*vap.d/(vap.d + aux.x*(liq.d - vap.d));
            aux.cv := Common.cv2Phase(liq, vap, aux.x, aux.T, p);
            aux.cp := liq.cp + aux.x*(vap.cp - liq.cp);
            aux.pt := liq.pt + aux.x*(vap.pt - liq.pt);
            aux.pd := liq.pd + aux.x*(vap.pd - liq.pd);
          elseif (aux.region == 5) then
            (aux.T, error) := BaseIF97.Inverses.tofph5(p = aux.p, h = aux.h, reldh = 1.0e-7);
            assert(error == 0, "Error in inverse iteration of steam tables");
            g := BaseIF97.Basic.g5(aux.p, aux.T);
            aux.s := aux.R_s*(g.tau*g.gtau - g.g);
            aux.rho := p/(aux.R_s*aux.T*g.pi*g.gpi);
            aux.vt := aux.R_s/p*(g.pi*g.gpi - g.tau*g.pi*g.gtaupi);
            aux.vp := aux.R_s*aux.T/(p*p)*g.pi*g.pi*g.gpipi;
            aux.pt := -g.p/g.T*(g.gpi - g.tau*g.gtaupi)/(g.gpipi*g.pi);
            aux.pd := -g.R_s*g.T*g.gpi*g.gpi/(g.gpipi);
            aux.cp := -aux.R_s*g.tau*g.tau*g.gtautau;
            aux.cv := aux.R_s*(-g.tau*g.tau*g.gtautau + ((g.gpi - g.tau*g.gtaupi)*(g.gpi - g.tau*g.gtaupi)/g.gpipi));
            aux.dpT := -aux.vt/aux.vp;
          else
            assert(false, "Error in region computation of IF97 steam tables" + "(p = " + String(p) + ", h = " + String(h) + ")");
          end if;
        end waterBaseProp_ph;

        function rho_props_ph
          input SI.Pressure p;
          input SI.SpecificEnthalpy h;
          input Common.IF97BaseTwoPhase properties;
          output SI.Density rho;
        algorithm
          rho := properties.rho;
          annotation(derivative(noDerivative = properties) = rho_ph_der, Inline = false, LateInline = true);
        end rho_props_ph;

        function rho_ph
          input SI.Pressure p;
          input SI.SpecificEnthalpy h;
          input Integer phase = 0;
          input Integer region = 0;
          output SI.Density rho;
        algorithm
          rho := rho_props_ph(p, h, waterBaseProp_ph(p, h, phase, region));
          annotation(Inline = true);
        end rho_ph;

        function rho_ph_der
          input SI.Pressure p;
          input SI.SpecificEnthalpy h;
          input Common.IF97BaseTwoPhase properties;
          input Real p_der;
          input Real h_der;
          output Real rho_der;
        algorithm
          if (properties.region == 4) then
            rho_der := (properties.rho*(properties.rho*properties.cv/properties.dpT + 1.0)/(properties.dpT*properties.T))*p_der + (-properties.rho*properties.rho/(properties.dpT*properties.T))*h_der;
          elseif (properties.region == 3) then
            rho_der := ((properties.rho*(properties.cv*properties.rho + properties.pt))/(properties.rho*properties.rho*properties.pd*properties.cv + properties.T*properties.pt*properties.pt))*p_der + (-properties.rho*properties.rho*properties.pt/(properties.rho*properties.rho*properties.pd*properties.cv + properties.T*properties.pt*properties.pt))*h_der;
          else
            rho_der := (-properties.rho*properties.rho*(properties.vp*properties.cp - properties.vt/properties.rho + properties.T*properties.vt*properties.vt)/properties.cp)*p_der + (-properties.rho*properties.rho*properties.vt/(properties.cp))*h_der;
          end if;
        end rho_ph_der;

        function T_props_ph
          input SI.Pressure p;
          input SI.SpecificEnthalpy h;
          input Common.IF97BaseTwoPhase properties;
          output SI.Temperature T;
        algorithm
          T := properties.T;
          annotation(derivative(noDerivative = properties) = T_ph_der, Inline = false, LateInline = true);
        end T_props_ph;

        function T_ph
          input SI.Pressure p;
          input SI.SpecificEnthalpy h;
          input Integer phase = 0;
          input Integer region = 0;
          output SI.Temperature T;
        algorithm
          T := T_props_ph(p, h, waterBaseProp_ph(p, h, phase, region));
          annotation(Inline = true);
        end T_ph;

        function T_ph_der
          input SI.Pressure p;
          input SI.SpecificEnthalpy h;
          input Common.IF97BaseTwoPhase properties;
          input Real p_der;
          input Real h_der;
          output Real T_der;
        algorithm
          if (properties.region == 4) then
            T_der := 1/properties.dpT*p_der;
          elseif (properties.region == 3) then
            T_der := ((-properties.rho*properties.pd + properties.T*properties.pt)/(properties.rho*properties.rho*properties.pd*properties.cv + properties.T*properties.pt*properties.pt))*p_der + ((properties.rho*properties.rho*properties.pd)/(properties.rho*properties.rho*properties.pd*properties.cv + properties.T*properties.pt*properties.pt))*h_der;
          else
            T_der := ((-1/properties.rho + properties.T*properties.vt)/properties.cp)*p_der + (1/properties.cp)*h_der;
          end if;
        end T_ph_der;

        function cp_props_ph
          input SI.Pressure p;
          input SI.SpecificEnthalpy h;
          input Common.IF97BaseTwoPhase aux;
          output SI.SpecificHeatCapacity cp;
        algorithm
          cp := aux.cp;
          annotation(Inline = false, LateInline = true);
        end cp_props_ph;

        function cp_ph
          input SI.Pressure p;
          input SI.SpecificEnthalpy h;
          input Integer phase = 0;
          input Integer region = 0;
          output SI.SpecificHeatCapacity cp;
        algorithm
          cp := cp_props_ph(p, h, waterBaseProp_ph(p, h, phase, region));
          annotation(Inline = true);
        end cp_ph;

        function waterBaseProp_pT
          input SI.Pressure p;
          input SI.Temperature T;
          input Integer region = 0;
          output Common.IF97BaseTwoPhase aux;
        protected
          Common.GibbsDerivs g;
          Common.HelmholtzDerivs f;
          Integer error;
        algorithm
          aux.phase := 1;
          aux.region := if region == 0 then BaseIF97.Regions.region_pT(p = p, T = T) else region;
          aux.R_s := BaseIF97.data.RH2O;
          aux.p := p;
          aux.T := T;
          aux.vt := 0.0;
          aux.vp := 0.0;
          if (aux.region == 1) then
            g := BaseIF97.Basic.g1(p, T);
            aux.h := aux.R_s*aux.T*g.tau*g.gtau;
            aux.s := aux.R_s*(g.tau*g.gtau - g.g);
            aux.rho := p/(aux.R_s*T*g.pi*g.gpi);
            aux.vt := aux.R_s/p*(g.pi*g.gpi - g.tau*g.pi*g.gtaupi);
            aux.vp := aux.R_s*T/(p*p)*g.pi*g.pi*g.gpipi;
            aux.cp := -aux.R_s*g.tau*g.tau*g.gtautau;
            aux.cv := aux.R_s*(-g.tau*g.tau*g.gtautau + ((g.gpi - g.tau*g.gtaupi)*(g.gpi - g.tau*g.gtaupi)/g.gpipi));
            aux.x := 0.0;
            aux.dpT := -aux.vt/aux.vp;
            aux.pt := -g.p/g.T*(g.gpi - g.tau*g.gtaupi)/(g.gpipi*g.pi);
            aux.pd := -g.R_s*g.T*g.gpi*g.gpi/(g.gpipi);
          elseif (aux.region == 2) then
            g := BaseIF97.Basic.g2(p, T);
            aux.h := aux.R_s*aux.T*g.tau*g.gtau;
            aux.s := aux.R_s*(g.tau*g.gtau - g.g);
            aux.rho := p/(aux.R_s*T*g.pi*g.gpi);
            aux.vt := aux.R_s/p*(g.pi*g.gpi - g.tau*g.pi*g.gtaupi);
            aux.vp := aux.R_s*T/(p*p)*g.pi*g.pi*g.gpipi;
            aux.pt := -g.p/g.T*(g.gpi - g.tau*g.gtaupi)/(g.gpipi*g.pi);
            aux.pd := -g.R_s*g.T*g.gpi*g.gpi/(g.gpipi);
            aux.cp := -aux.R_s*g.tau*g.tau*g.gtautau;
            aux.cv := aux.R_s*(-g.tau*g.tau*g.gtautau + ((g.gpi - g.tau*g.gtaupi)*(g.gpi - g.tau*g.gtaupi)/g.gpipi));
            aux.x := 1.0;
            aux.dpT := -aux.vt/aux.vp;
          elseif (aux.region == 3) then
            (aux.rho, error) := BaseIF97.Inverses.dofpt3(p = p, T = T, delp = 1.0e-7);
            f := BaseIF97.Basic.f3(aux.rho, T);
            aux.h := aux.R_s*T*(f.tau*f.ftau + f.delta*f.fdelta);
            aux.s := aux.R_s*(f.tau*f.ftau - f.f);
            aux.pd := aux.R_s*T*f.delta*(2.0*f.fdelta + f.delta*f.fdeltadelta);
            aux.pt := aux.R_s*aux.rho*f.delta*(f.fdelta - f.tau*f.fdeltatau);
            aux.cv := aux.R_s*(-f.tau*f.tau*f.ftautau);
            aux.x := 0.0;
            aux.dpT := aux.pt;
          elseif (aux.region == 5) then
            g := BaseIF97.Basic.g5(p, T);
            aux.h := aux.R_s*aux.T*g.tau*g.gtau;
            aux.s := aux.R_s*(g.tau*g.gtau - g.g);
            aux.rho := p/(aux.R_s*T*g.pi*g.gpi);
            aux.vt := aux.R_s/p*(g.pi*g.gpi - g.tau*g.pi*g.gtaupi);
            aux.vp := aux.R_s*T/(p*p)*g.pi*g.pi*g.gpipi;
            aux.pt := -g.p/g.T*(g.gpi - g.tau*g.gtaupi)/(g.gpipi*g.pi);
            aux.pd := -g.R_s*g.T*g.gpi*g.gpi/(g.gpipi);
            aux.cp := -aux.R_s*g.tau*g.tau*g.gtautau;
            aux.cv := aux.R_s*(-g.tau*g.tau*g.gtautau + ((g.gpi - g.tau*g.gtaupi)*(g.gpi - g.tau*g.gtaupi)/g.gpipi));
            aux.x := 1.0;
            aux.dpT := -aux.vt/aux.vp;
          else
            assert(false, "Error in region computation of IF97 steam tables" + "(p = " + String(p) + ", T = " + String(T) + ")");
          end if;
        end waterBaseProp_pT;

        function rho_props_pT
          input SI.Pressure p;
          input SI.Temperature T;
          input Common.IF97BaseTwoPhase aux;
          output SI.Density rho;
        algorithm
          rho := aux.rho;
          annotation(derivative(noDerivative = aux) = rho_pT_der, Inline = false, LateInline = true);
        end rho_props_pT;

        function rho_pT
          input SI.Pressure p;
          input SI.Temperature T;
          input Integer region = 0;
          output SI.Density rho;
        algorithm
          rho := rho_props_pT(p, T, waterBaseProp_pT(p, T, region));
          annotation(Inline = true);
        end rho_pT;

        function h_props_pT
          input SI.Pressure p;
          input SI.Temperature T;
          input Common.IF97BaseTwoPhase aux;
          output SI.SpecificEnthalpy h;
        algorithm
          h := aux.h;
          annotation(derivative(noDerivative = aux) = h_pT_der, Inline = false, LateInline = true);
        end h_props_pT;

        function h_pT
          input SI.Pressure p;
          input SI.Temperature T;
          input Integer region = 0;
          output SI.SpecificEnthalpy h;
        algorithm
          h := h_props_pT(p, T, waterBaseProp_pT(p, T, region));
          annotation(Inline = true);
        end h_pT;

        function h_pT_der
          input SI.Pressure p;
          input SI.Temperature T;
          input Common.IF97BaseTwoPhase aux;
          input Real p_der;
          input Real T_der;
          output Real h_der;
        algorithm
          if (aux.region == 3) then
            h_der := ((-aux.rho*aux.pd + T*aux.pt)/(aux.rho*aux.rho*aux.pd))*p_der + ((aux.rho*aux.rho*aux.pd*aux.cv + aux.T*aux.pt*aux.pt)/(aux.rho*aux.rho*aux.pd))*T_der;
          else
            h_der := (1/aux.rho - aux.T*aux.vt)*p_der + aux.cp*T_der;
          end if;
        end h_pT_der;

        function rho_pT_der
          input SI.Pressure p;
          input SI.Temperature T;
          input Common.IF97BaseTwoPhase aux;
          input Real p_der;
          input Real T_der;
          output Real rho_der;
        algorithm
          if (aux.region == 3) then
            rho_der := (1/aux.pd)*p_der - (aux.pt/aux.pd)*T_der;
          else
            rho_der := (-aux.rho*aux.rho*aux.vp)*p_der + (-aux.rho*aux.rho*aux.vt)*T_der;
          end if;
        end rho_pT_der;

        function cp_props_pT
          input SI.Pressure p;
          input SI.Temperature T;
          input Common.IF97BaseTwoPhase aux;
          output SI.SpecificHeatCapacity cp;
        algorithm
          cp := if aux.region == 3 then (aux.rho*aux.rho*aux.pd*aux.cv + aux.T*aux.pt*aux.pt)/(aux.rho*aux.rho*aux.pd) else aux.cp;
          annotation(Inline = false, LateInline = true);
        end cp_props_pT;

        function cp_pT
          input SI.Pressure p;
          input SI.Temperature T;
          input Integer region = 0;
          output SI.SpecificHeatCapacity cp;
        algorithm
          cp := cp_props_pT(p, T, waterBaseProp_pT(p, T, region));
          annotation(Inline = true);
        end cp_pT;

        function waterBaseProp_dT
          input SI.Density rho;
          input SI.Temperature T;
          input Integer phase = 0;
          input Integer region = 0;
          output Common.IF97BaseTwoPhase aux;
        protected
          SI.SpecificEnthalpy h_liq;
          SI.Density d_liq;
          SI.SpecificEnthalpy h_vap;
          SI.Density d_vap;
          Common.GibbsDerivs g;
          Common.HelmholtzDerivs f;
          Modelica.Media.Common.PhaseBoundaryProperties liq;
          Modelica.Media.Common.PhaseBoundaryProperties vap;
          Modelica.Media.Common.GibbsDerivs gl;
          Modelica.Media.Common.GibbsDerivs gv;
          Modelica.Media.Common.HelmholtzDerivs fl;
          Modelica.Media.Common.HelmholtzDerivs fv;
          Integer error;
        algorithm
          aux.region := if region == 0 then (if phase == 2 then 4 else BaseIF97.Regions.region_dT(d = rho, T = T, phase = phase)) else region;
          aux.phase := if aux.region == 4 then 2 else 1;
          aux.R_s := BaseIF97.data.RH2O;
          aux.rho := rho;
          aux.T := T;
          aux.vt := 0.0;
          aux.vp := 0.0;
          if (aux.region == 1) then
            (aux.p, error) := BaseIF97.Inverses.pofdt125(d = rho, T = T, reldd = 1.0e-8, region = 1);
            g := BaseIF97.Basic.g1(aux.p, T);
            aux.h := aux.R_s*aux.T*g.tau*g.gtau;
            aux.s := aux.R_s*(g.tau*g.gtau - g.g);
            aux.rho := aux.p/(aux.R_s*T*g.pi*g.gpi);
            aux.vt := aux.R_s/aux.p*(g.pi*g.gpi - g.tau*g.pi*g.gtaupi);
            aux.vp := aux.R_s*T/(aux.p*aux.p)*g.pi*g.pi*g.gpipi;
            aux.pt := -g.p/g.T*(g.gpi - g.tau*g.gtaupi)/(g.gpipi*g.pi);
            aux.pd := -g.R_s*g.T*g.gpi*g.gpi/(g.gpipi);
            aux.cp := -aux.R_s*g.tau*g.tau*g.gtautau;
            aux.cv := aux.R_s*(-g.tau*g.tau*g.gtautau + ((g.gpi - g.tau*g.gtaupi)*(g.gpi - g.tau*g.gtaupi)/g.gpipi));
            aux.x := 0.0;
          elseif (aux.region == 2) then
            (aux.p, error) := BaseIF97.Inverses.pofdt125(d = rho, T = T, reldd = 1.0e-8, region = 2);
            g := BaseIF97.Basic.g2(aux.p, T);
            aux.h := aux.R_s*aux.T*g.tau*g.gtau;
            aux.s := aux.R_s*(g.tau*g.gtau - g.g);
            aux.rho := aux.p/(aux.R_s*T*g.pi*g.gpi);
            aux.vt := aux.R_s/aux.p*(g.pi*g.gpi - g.tau*g.pi*g.gtaupi);
            aux.vp := aux.R_s*T/(aux.p*aux.p)*g.pi*g.pi*g.gpipi;
            aux.pt := -g.p/g.T*(g.gpi - g.tau*g.gtaupi)/(g.gpipi*g.pi);
            aux.pd := -g.R_s*g.T*g.gpi*g.gpi/(g.gpipi);
            aux.cp := -aux.R_s*g.tau*g.tau*g.gtautau;
            aux.cv := aux.R_s*(-g.tau*g.tau*g.gtautau + ((g.gpi - g.tau*g.gtaupi)*(g.gpi - g.tau*g.gtaupi)/g.gpipi));
            aux.x := 1.0;
          elseif (aux.region == 3) then
            f := BaseIF97.Basic.f3(rho, T);
            aux.p := aux.R_s*rho*T*f.delta*f.fdelta;
            aux.h := aux.R_s*T*(f.tau*f.ftau + f.delta*f.fdelta);
            aux.s := aux.R_s*(f.tau*f.ftau - f.f);
            aux.pd := aux.R_s*T*f.delta*(2.0*f.fdelta + f.delta*f.fdeltadelta);
            aux.pt := aux.R_s*rho*f.delta*(f.fdelta - f.tau*f.fdeltatau);
            aux.cv := aux.R_s*(-f.tau*f.tau*f.ftautau);
            aux.cp := (aux.rho*aux.rho*aux.pd*aux.cv + aux.T*aux.pt*aux.pt)/(aux.rho*aux.rho*aux.pd);
            aux.x := 0.0;
          elseif (aux.region == 4) then
            aux.p := BaseIF97.Basic.psat(T);
            d_liq := rhol_T(T);
            d_vap := rhov_T(T);
            h_liq := hl_p(aux.p);
            h_vap := hv_p(aux.p);
            aux.x := if (d_vap <> d_liq) then (1/rho - 1/d_liq)/(1/d_vap - 1/d_liq) else 1.0;
            aux.h := h_liq + aux.x*(h_vap - h_liq);
            if T < BaseIF97.data.TLIMIT1 then
              gl := BaseIF97.Basic.g1(aux.p, T);
              gv := BaseIF97.Basic.g2(aux.p, T);
              liq := Common.gibbsToBoundaryProps(gl);
              vap := Common.gibbsToBoundaryProps(gv);
            else
              fl := BaseIF97.Basic.f3(d_liq, T);
              fv := BaseIF97.Basic.f3(d_vap, T);
              liq := Common.helmholtzToBoundaryProps(fl);
              vap := Common.helmholtzToBoundaryProps(fv);
            end if;
            aux.dpT := if (liq.d <> vap.d) then (vap.s - liq.s)*liq.d*vap.d/(liq.d - vap.d) else BaseIF97.Basic.dptofT(aux.T);
            aux.s := liq.s + aux.x*(vap.s - liq.s);
            aux.cv := Common.cv2Phase(liq, vap, aux.x, aux.T, aux.p);
            aux.cp := liq.cp + aux.x*(vap.cp - liq.cp);
            aux.pt := liq.pt + aux.x*(vap.pt - liq.pt);
            aux.pd := liq.pd + aux.x*(vap.pd - liq.pd);
          elseif (aux.region == 5) then
            (aux.p, error) := BaseIF97.Inverses.pofdt125(d = rho, T = T, reldd = 1.0e-8, region = 5);
            g := BaseIF97.Basic.g2(aux.p, T);
            aux.h := aux.R_s*aux.T*g.tau*g.gtau;
            aux.s := aux.R_s*(g.tau*g.gtau - g.g);
            aux.rho := aux.p/(aux.R_s*T*g.pi*g.gpi);
            aux.vt := aux.R_s/aux.p*(g.pi*g.gpi - g.tau*g.pi*g.gtaupi);
            aux.vp := aux.R_s*T/(aux.p*aux.p)*g.pi*g.pi*g.gpipi;
            aux.pt := -g.p/g.T*(g.gpi - g.tau*g.gtaupi)/(g.gpipi*g.pi);
            aux.pd := -g.R_s*g.T*g.gpi*g.gpi/(g.gpipi);
            aux.cp := -aux.R_s*g.tau*g.tau*g.gtautau;
            aux.cv := aux.R_s*(-g.tau*g.tau*g.gtautau + ((g.gpi - g.tau*g.gtaupi)*(g.gpi - g.tau*g.gtaupi)/g.gpipi));
          else
            assert(false, "Error in region computation of IF97 steam tables" + "(rho = " + String(rho) + ", T = " + String(T) + ")");
          end if;
        end waterBaseProp_dT;

        function h_props_dT
          input SI.Density d;
          input SI.Temperature T;
          input Common.IF97BaseTwoPhase aux;
          output SI.SpecificEnthalpy h;
        algorithm
          h := aux.h;
          annotation(derivative(noDerivative = aux) = h_dT_der, Inline = false, LateInline = true);
        end h_props_dT;

        function h_dT
          input SI.Density d;
          input SI.Temperature T;
          input Integer phase = 0;
          input Integer region = 0;
          output SI.SpecificEnthalpy h;
        algorithm
          h := h_props_dT(d, T, waterBaseProp_dT(d, T, phase, region));
          annotation(Inline = true);
        end h_dT;

        function h_dT_der
          input SI.Density d;
          input SI.Temperature T;
          input Common.IF97BaseTwoPhase aux;
          input Real d_der;
          input Real T_der;
          output Real h_der;
        algorithm
          if (aux.region == 3) then
            h_der := ((-d*aux.pd + T*aux.pt)/(d*d))*d_der + ((aux.cv*d + aux.pt)/d)*T_der;
          elseif (aux.region == 4) then
            h_der := T*aux.dpT/(d*d)*d_der + ((aux.cv*d + aux.dpT)/d)*T_der;
          else
            h_der := (-(-1/d + T*aux.vt)/(d*d*aux.vp))*d_der + ((aux.vp*aux.cp - aux.vt/d + T*aux.vt*aux.vt)/aux.vp)*T_der;
          end if;
        end h_dT_der;

        function p_props_dT
          input SI.Density d;
          input SI.Temperature T;
          input Common.IF97BaseTwoPhase aux;
          output SI.Pressure p;
        algorithm
          p := aux.p;
          annotation(derivative(noDerivative = aux) = p_dT_der, Inline = false, LateInline = true);
        end p_props_dT;

        function p_dT
          input SI.Density d;
          input SI.Temperature T;
          input Integer phase = 0;
          input Integer region = 0;
          output SI.Pressure p;
        algorithm
          p := p_props_dT(d, T, waterBaseProp_dT(d, T, phase, region));
          annotation(Inline = true);
        end p_dT;

        function p_dT_der
          input SI.Density d;
          input SI.Temperature T;
          input Common.IF97BaseTwoPhase aux;
          input Real d_der;
          input Real T_der;
          output Real p_der;
        algorithm
          if (aux.region == 3) then
            p_der := aux.pd*d_der + aux.pt*T_der;
          elseif (aux.region == 4) then
            p_der := aux.dpT*T_der;
          else
            p_der := (-1/(d*d*aux.vp))*d_der + (-aux.vt/aux.vp)*T_der;
          end if;
        end p_dT_der;

        function cp_props_dT
          input SI.Density d;
          input SI.Temperature T;
          input Common.IF97BaseTwoPhase aux;
          output SI.SpecificHeatCapacity cp;
        algorithm
          cp := aux.cp;
          annotation(Inline = false, LateInline = true);
        end cp_props_dT;

        function cp_dT
          input SI.Density d;
          input SI.Temperature T;
          input Integer phase = 0;
          input Integer region = 0;
          output SI.SpecificHeatCapacity cp;
        algorithm
          cp := cp_props_dT(d, T, waterBaseProp_dT(d, T, phase, region));
          annotation(Inline = true);
        end cp_dT;

        function hl_p = BaseIF97.Regions.hl_p;
        function hv_p = BaseIF97.Regions.hv_p;
        function rhol_T = BaseIF97.Regions.rhol_T;
        function rhov_T = BaseIF97.Regions.rhov_T;
        function rhol_p = BaseIF97.Regions.rhol_p;
        function rhov_p = BaseIF97.Regions.rhov_p;
      end IF97_Utilities;
    end Water;
  end Media;

  package Math
    package Icons
      partial function AxisLeft end AxisLeft;

      partial function AxisCenter end AxisCenter;
    end Icons;

    function asin
      input Real u;
      output Modelica.Units.SI.Angle y;
      external "builtin" y = asin(u);
    end asin;

    function acos
      input Real u;
      output Modelica.Units.SI.Angle y;
      external "builtin" y = acos(u);
    end acos;

    function exp
      input Real u;
      output Real y;
      external "builtin" y = exp(u);
    end exp;

    function log
      input Real u;
      output Real y;
      external "builtin" y = log(u);
    end log;
  end Math;

  package Utilities
    package Types end Types;

    package Internal
      import Modelica.Units.SI;
    end Internal;
  end Utilities;

  package Constants
    import Modelica.Units.SI;
    import Modelica.Units.NonSI;
    final constant Real e = Modelica.Math.exp(1.0);
    final constant Real pi = 2*Modelica.Math.asin(1.0);
    final constant Real D2R = pi/180;
    final constant Real R2D = 180/pi;
    final constant Real gamma = 0.57721566490153286061;
    final constant Real eps = ModelicaServices.Machine.eps;
    final constant Real small = ModelicaServices.Machine.small;
    final constant Real inf = ModelicaServices.Machine.inf;
    final constant Integer Integer_inf = ModelicaServices.Machine.Integer_inf;
    final constant SI.Velocity c = 299792458;
    final constant SI.Acceleration g_n = 9.80665;
    final constant Real G(final unit = "m3/(kg.s2)") = 6.67430e-11;
    final constant SI.ElectricCharge q = 1.602176634e-19;
    final constant SI.FaradayConstant F = q*N_A;
    final constant Real h(final unit = "J.s") = 6.62607015e-34;
    final constant Real k(final unit = "J/K") = 1.380649e-23;
    final constant Real R(final unit = "J/(mol.K)") = k*N_A;
    final constant Real sigma(final unit = "W/(m2.K4)") = 2*pi^5*k^4/(15*h^3*c^2);
    final constant Real N_A(final unit = "1/mol") = 6.02214076e23;
    final constant Real mu_0(final unit = "N/A2") = 4*pi*1.00000000055e-7;
    final constant Real epsilon_0(final unit = "F/m") = 1/(mu_0*c*c);
    final constant NonSI.Temperature_degC T_zero = -273.15;
  end Constants;

  package Icons
    partial package Package end Package;

    partial package VariantsPackage end VariantsPackage;

    partial package InterfacesPackage end InterfacesPackage;

    partial package UtilitiesPackage end UtilitiesPackage;

    partial package TypesPackage end TypesPackage;

    partial package FunctionsPackage end FunctionsPackage;

    partial package IconsPackage end IconsPackage;

    partial package InternalPackage end InternalPackage;

    partial package MaterialPropertiesPackage end MaterialPropertiesPackage;

    partial function Function end Function;

    partial record Record end Record;
  end Icons;

  package Units
    package SI
      type Angle = Real(final quantity = "Angle", final unit = "rad", displayUnit = "deg");
      type Area = Real(final quantity = "Area", final unit = "m2");
      type Volume = Real(final quantity = "Volume", final unit = "m3");
      type Velocity = Real(final quantity = "Velocity", final unit = "m/s");
      type Acceleration = Real(final quantity = "Acceleration", final unit = "m/s2");
      type Mass = Real(quantity = "Mass", final unit = "kg", min = 0);
      type Density = Real(final quantity = "Density", final unit = "kg/m3", displayUnit = "g/cm3", min = 0.0);
      type SpecificVolume = Real(final quantity = "SpecificVolume", final unit = "m3/kg", min = 0.0);
      type Pressure = Real(final quantity = "Pressure", final unit = "Pa", displayUnit = "bar");
      type AbsolutePressure = Pressure(min = 0.0, nominal = 1e5);
      type DynamicViscosity = Real(final quantity = "DynamicViscosity", final unit = "Pa.s", min = 0);
      type Energy = Real(final quantity = "Energy", final unit = "J");
      type Power = Real(final quantity = "Power", final unit = "W");
      type MassFlowRate = Real(quantity = "MassFlowRate", final unit = "kg/s");
      type MomentumFlux = Real(final quantity = "MomentumFlux", final unit = "N");
      type ThermodynamicTemperature = Real(final quantity = "ThermodynamicTemperature", final unit = "K", min = 0.0, start = 288.15, nominal = 300, displayUnit = "degC") annotation(absoluteValue = true);
      type Temperature = ThermodynamicTemperature;
      type Heat = Real(final quantity = "Energy", final unit = "J");
      type ThermalConductivity = Real(final quantity = "ThermalConductivity", final unit = "W/(m.K)");
      type SpecificHeatCapacity = Real(final quantity = "SpecificHeatCapacity", final unit = "J/(kg.K)");
      type RatioOfSpecificHeatCapacities = Real(final quantity = "RatioOfSpecificHeatCapacities", final unit = "1");
      type Entropy = Real(final quantity = "Entropy", final unit = "J/K");
      type SpecificEntropy = Real(final quantity = "SpecificEntropy", final unit = "J/(kg.K)");
      type SpecificEnergy = Real(final quantity = "SpecificEnergy", final unit = "J/kg");
      type SpecificInternalEnergy = SpecificEnergy;
      type SpecificEnthalpy = SpecificEnergy;
      type DerPressureByDensity = Real(final unit = "Pa.m3/kg");
      type DerPressureByTemperature = Real(final unit = "Pa/K");
      type ElectricCharge = Real(final quantity = "ElectricCharge", final unit = "C");
      type VelocityOfSound = Real(final quantity = "Velocity", final unit = "m/s");
      type AmountOfSubstance = Real(final quantity = "AmountOfSubstance", final unit = "mol", min = 0);
      type MolarMass = Real(final quantity = "MolarMass", final unit = "kg/mol", min = 0);
      type MolarVolume = Real(final quantity = "MolarVolume", final unit = "m3/mol", min = 0);
      type MassFraction = Real(final quantity = "MassFraction", final unit = "1", min = 0, max = 1);
      type MoleFraction = Real(final quantity = "MoleFraction", final unit = "1", min = 0, max = 1);
      type FaradayConstant = Real(final quantity = "FaradayConstant", final unit = "C/mol");
    end SI;

    package NonSI
      type Temperature_degC = Real(final quantity = "ThermodynamicTemperature", final unit = "degC") annotation(absoluteValue = true);
      type Pressure_bar = Real(final quantity = "Pressure", final unit = "bar");
    end NonSI;

    package Conversions
      function to_degC
        input SI.Temperature Kelvin;
        output Modelica.Units.NonSI.Temperature_degC Celsius;
      algorithm
        Celsius := Kelvin + Modelica.Constants.T_zero;
        annotation(Inline = true);
      end to_degC;

      function from_degC
        input Modelica.Units.NonSI.Temperature_degC Celsius;
        output SI.Temperature Kelvin;
      algorithm
        Kelvin := Celsius - Modelica.Constants.T_zero;
        annotation(Inline = true);
      end from_degC;

      function to_bar
        input SI.Pressure Pa;
        output Modelica.Units.NonSI.Pressure_bar bar;
      algorithm
        bar := Pa/1e5;
        annotation(Inline = true);
      end to_bar;
    end Conversions;

    package Icons
      partial function Conversion end Conversion;
    end Icons;
  end Units;
  annotation(version = "4.0.0", versionDate = "2020-06-04", dateModified = "2020-06-04 11:00:00Z");
end Modelica;

model plant
  import ThermoS.Uops.*;
  import ThermoS.Types.*;
  import Water = Modelica.Media.Water.StandardWater;
  import ThermoS.Uops.Tanks.SimpleTank;
  SimpleTank tank1(redeclare package Medium = Water);
  SimpleTank tank2(redeclare package Medium = Water);
  Feed feed1(redeclare package Medium = Water);
  Feed feed2(redeclare package Medium = Water);
  Product prod1(redeclare package Medium = Water);
  Product prod2(redeclare package Medium = Water);
  Real U;
  Real T1;
  Real T2;
initial algorithm
  tank1.mf := 50;
  tank1.Tf := 20 + 273;
  tank2.mf := 100;
  tank2.Tf := 10 + 273;
equation
  connect(feed1.outlet, tank1.inlet);
  connect(tank1.outlet, prod1.inlet);
  connect(feed2.outlet, tank2.inlet);
  connect(tank2.outlet, prod2.inlet);
  feed1.mdot = 2.0/60;
  feed2.mdot = 2.0/60;
  feed1.T = tank2.Tf;
  feed1.Xi = tank2.outlet.Xi_outflow;
  feed2.T = tank1.Tf;
  feed2.Xi = tank1.outlet.Xi_outflow;
  tank1.Q_in = U*0.5*((100 + 273) - tank1.Tf);
  tank2.Q_in = U*0.5*((120 + 273) - tank2.Tf);
  tank1.tvol = 1000;
  tank2.tvol = 1000;
  tank1.inlet.p = tank1.Pa;
  tank2.inlet.p = tank2.Pa;
  prod1.inlet.m_flow = feed1.mdot;
  prod2.inlet.m_flow = feed2.mdot;
  prod1.T = 300;
  prod1.p = 1e5;
  prod1.Xi = fill(0.0, 0);
  prod2.T = 300;
  prod2.p = 1e5;
  prod2.Xi = fill(0.0, 0);
  U = (4000*4.17)/60;
  tank1.Pa = 1e5;
  tank2.Pa = 1e5;
  T1 = tank1.Tf - 273;
  T2 = tank2.Tf - 273;
end plant;

package ThermoS
  package Uops
    package Tanks
      model SimpleTank
        replaceable package Medium = PartialMixtureMedium;
        FluidPort inlet(redeclare package Medium = Medium);
        FluidPort outlet(redeclare package Medium = Medium);
        Heat Q_in;
        Temperature Tf;
        Temperature Pa;
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
    end Tanks;

    model Product
      replaceable package Medium = PartialMixtureMedium;
      FluidPort inlet(redeclare package Medium = Medium);
      Medium.AbsolutePressure p;
      Medium.Temperature T;
      Medium.MassFraction[Medium.nXi] Xi;
      Medium.BaseProperties medium;
    equation
      medium.p = inlet.p;
      medium.T = T;
      medium.Xi = Xi;
      inlet.h_outflow = medium.h;
      inlet.Xi_outflow = Xi;
    end Product;

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
    end Feed;

    import Modelica.Media.Interfaces.PartialMixtureMedium;
    import Modelica.Fluid.Interfaces.FluidPort;
    import Modelica.Fluid.Utilities.*;
    import Modelica.Units.SI.*;
    import ThermoS.Types.*;
  end Uops;

  package Types end Types;

  package Media end Media;

  package Math end Math;
end ThermoS;

model plant_total
  extends plant;
end plant_total;
