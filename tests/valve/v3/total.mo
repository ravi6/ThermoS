model plant  
  package Gas = ThermoS.Media.MyGas(MassFlowRate(min = -10, max = 20));
  constant Real[2] AirComp = {0.767, 0.233};
  .ThermoS.Uops.Reservoir src(redeclare package Medium = Gas, p = 2.0e5, each T = 300, Xi = AirComp);
  .ThermoS.Uops.Reservoir sink(redeclare each package Medium = Gas, p = 2.0e5, T = 300, Xi = AirComp);
  .ThermoS.Uops.Valves.Valve v1(redeclare package Medium = Gas, cv = 0.004 / sqrt(0.5e5), Compressible = true);
  .ThermoS.Uops.Valves.Valve v2(redeclare package Medium = Gas, cv = 0.004 / sqrt(0.5e5), Compressible = true);
equation
  v1.po = 100;
  v2.po = 0;
  connect(src.port, v1.inlet);
  connect(src.port, v2.inlet);
  connect(v1.outlet, sink.port);
  connect(v2.outlet, sink.port);
end plant;

package ThermoS  "A Modelica Package for Process Simulations" 
  package Uops  "Unit Operations in ThermoS Package" 
    package Valves  "Contains Valves and its interfaces " 
      partial model partialValve  
        replaceable package Medium = .Modelica.Media.Interfaces.PartialMixtureMedium;
        .Modelica.Fluid.Interfaces.FluidPort inlet(redeclare package Medium = Medium);
        .Modelica.Fluid.Interfaces.FluidPort outlet(redeclare package Medium = Medium);
        parameter Real cv = 1.0 / sqrt(1e5);
        parameter Real dpTol = 100;
        Medium.BaseProperties med(preferredMediumStates = true);
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
        parameter .ThermoS.Types.Fraction pratChoke = 0.5;
        parameter Boolean Compressible = true;
        .ThermoS.Types.Percent po(start = 1.0);
        .ThermoS.Types.Fraction charF(start = 1.0);
        .ThermoS.Types.Fraction prat(start = 1.0);
      equation
        if vchar == Vchar.Linear then
          charF = po / 100;
        elseif vchar == Vchar.FastActing then
          charF = (po / 100) ^ 0.5;
        elseif vchar == Vchar.EquiPercent then
          charF = 35 ^ (po / 100 - 1);
        end if;
        prat = min(inlet.p, outlet.p) / max(inlet.p, outlet.p);
        if Compressible then
          inlet.m_flow = cv * charF * sqrt(abs(med.d)) * sign(inlet.p - outlet.p) * sqrt(max(inlet.p, outlet.p)) * (if prat > 0.5 then .Modelica.Fluid.Utilities.regRoot(1 - prat, dpTol) else sqrt(0.5));
        else
          inlet.m_flow = cv * charF * sqrt(med.d * inlet.p) * .Modelica.Fluid.Utilities.regRoot(1 - outlet.p / inlet.p, dpTol);
        end if;
      end Valve;

      type Vchar = enumeration(Linear "Linear Valve", FastActing "Fast Acting Valve", EquiPercent "Equi-Percent Valve") "Enumeration Defining Valve Behaviour";
    end Valves;

    model Reservoir  
      replaceable package Medium = .Modelica.Media.Interfaces.PartialMixtureMedium;
      .Modelica.Fluid.Interfaces.FluidPort port(redeclare package Medium = Medium);
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
  end Uops;

  package Types  "Contains all ThermoS Typed Variables" 
    type Percent = Real(unit = "%", min = 0, max = 100);
    type Fraction = Real(min = 0, max = 1.0);
  end Types;

  package Media  "All fluids that ThermoS knows about" 
    package MyGas  "Specific ideal Gas Mixture" 
      extends .Modelica.Media.IdealGases.Common.MixtureGasNasa(preferredMediumStates = true, data = {.Modelica.Media.IdealGases.Common.SingleGasesData.N2, .Modelica.Media.IdealGases.Common.SingleGasesData.O2, .Modelica.Media.IdealGases.Common.SingleGasesData.CO2}, fluidConstants = {.Modelica.Media.IdealGases.Common.FluidData.N2, .Modelica.Media.IdealGases.Common.FluidData.O2, .Modelica.Media.IdealGases.Common.FluidData.CO2}, substanceNames = {"Nitrogen", "Oxygen", "CarbonDioxide"}, reducedX = true, reference_X = {0.7, 0.2, 0.1}, Density(start = 1, nominal = 1), AbsolutePressure(start = 1e5, min = 1e3, max = 50e5, nominal = 1e5), Temperature(start = 300, min = 200, max = 2000, nominal = 300), ThermodynamicState(p(start = 1e5), T(start = 300), X(start = reference_X)), MassFraction(start = 0.333333));
    end MyGas;
  end Media;
end ThermoS;

package ModelicaServices  "ModelicaServices (OpenModelica implementation) - Models and functions used in the Modelica Standard Library requiring a tool specific implementation" 
  extends Modelica.Icons.Package;

  package Machine  
    extends Modelica.Icons.Package;
    final constant Real eps = 1.e-15 "Biggest number such that 1.0 + eps = 1.0";
    final constant Real small = 1.e-60 "Smallest number such that small and -small are representable on the machine";
    final constant Real inf = 1.e+60 "Biggest Real number such that inf and -inf are representable on the machine";
    final constant Integer Integer_inf = OpenModelica.Internal.Architecture.integerMax() "Biggest Integer number such that Integer_inf and -Integer_inf are representable on the machine";
  end Machine;
  annotation(Protection(access = Access.hide), version = "3.2.2", versionBuild = 0, versionDate = "2016-01-15", dateModified = "2016-01-15 08:44:41Z"); 
end ModelicaServices;

package Modelica  "Modelica Standard Library - Version 3.2.2" 
  extends Modelica.Icons.Package;

  package Fluid  "Library of 1-dim. thermo-fluid flow models using the Modelica.Media media description" 
    extends Modelica.Icons.Package;

    package Interfaces  "Interfaces for steady state and unsteady, mixed-phase, multi-substance, incompressible and compressible flow" 
      extends Modelica.Icons.InterfacesPackage;

      connector FluidPort  "Interface for quasi one-dimensional fluid flow in a piping network (incompressible or compressible, one or more phases, one or more substances)" 
        replaceable package Medium = Modelica.Media.Interfaces.PartialMedium "Medium model" annotation(choicesAllMatching = true);
        flow Medium.MassFlowRate m_flow "Mass flow rate from the connection point into the component";
        Medium.AbsolutePressure p "Thermodynamic pressure in the connection point";
        stream Medium.SpecificEnthalpy h_outflow "Specific thermodynamic enthalpy close to the connection point if m_flow < 0";
        stream Medium.MassFraction[Medium.nXi] Xi_outflow "Independent mixture mass fractions m_i/m close to the connection point if m_flow < 0";
        stream Medium.ExtraProperty[Medium.nC] C_outflow "Properties c_i/m close to the connection point if m_flow < 0";
      end FluidPort;
    end Interfaces;

    package Utilities  "Utility models to construct fluid components (should not be used directly)" 
      extends Modelica.Icons.UtilitiesPackage;

      function regRoot  "Anti-symmetric square root approximation with finite derivative in the origin" 
        extends Modelica.Icons.Function;
        input Real x;
        input Real delta = 0.01 "Range of significant deviation from sqrt(abs(x))*sgn(x)";
        output Real y;
      algorithm
        y := x / (x * x + delta * delta) ^ 0.25;
        annotation(derivative(zeroDerivative = delta) = regRoot_der); 
      end regRoot;

      function regRoot_der  "Derivative of regRoot" 
        extends Modelica.Icons.Function;
        input Real x;
        input Real delta = 0.01 "Range of significant deviation from sqrt(x)";
        input Real dx "Derivative of x";
        output Real dy;
      algorithm
        dy := dx * 0.5 * (x * x + 2 * delta * delta) / (x * x + delta * delta) ^ 1.25;
      end regRoot_der;
    end Utilities;
  end Fluid;

  package Media  "Library of media property models" 
    extends Modelica.Icons.Package;

    package Interfaces  "Interfaces for media models" 
      extends Modelica.Icons.InterfacesPackage;

      partial package PartialMedium  "Partial medium properties (base package of all media packages)" 
        extends Modelica.Media.Interfaces.Types;
        extends Modelica.Icons.MaterialPropertiesPackage;
        constant Modelica.Media.Interfaces.Choices.IndependentVariables ThermoStates "Enumeration type for independent variables";
        constant String mediumName = "unusablePartialMedium" "Name of the medium";
        constant String[:] substanceNames = {mediumName} "Names of the mixture substances. Set substanceNames={mediumName} if only one substance.";
        constant String[:] extraPropertiesNames = fill("", 0) "Names of the additional (extra) transported properties. Set extraPropertiesNames=fill(\"\",0) if unused";
        constant Boolean singleState "= true, if u and d are not a function of pressure";
        constant Boolean reducedX = true "= true if medium contains the equation sum(X) = 1.0; set reducedX=true if only one substance (see docu for details)";
        constant Boolean fixedX = false "= true if medium contains the equation X = reference_X";
        constant AbsolutePressure reference_p = 101325 "Reference pressure of Medium: default 1 atmosphere";
        constant MassFraction[nX] reference_X = fill(1 / nX, nX) "Default mass fractions of medium";
        constant AbsolutePressure p_default = 101325 "Default value for pressure of medium (for initialization)";
        constant Temperature T_default = Modelica.SIunits.Conversions.from_degC(20) "Default value for temperature of medium (for initialization)";
        constant MassFraction[nX] X_default = reference_X "Default value for mass fractions of medium (for initialization)";
        final constant Integer nS = size(substanceNames, 1) "Number of substances" annotation(Evaluate = true);
        constant Integer nX = nS "Number of mass fractions" annotation(Evaluate = true);
        constant Integer nXi = if fixedX then 0 else if reducedX then nS - 1 else nS "Number of structurally independent mass fractions (see docu for details)" annotation(Evaluate = true);
        final constant Integer nC = size(extraPropertiesNames, 1) "Number of extra (outside of standard mass-balance) transported properties" annotation(Evaluate = true);
        replaceable record FluidConstants = Modelica.Media.Interfaces.Types.Basic.FluidConstants "Critical, triple, molecular and other standard data of fluid";

        replaceable record ThermodynamicState  "Minimal variable set that is available as input argument to every medium function" 
          extends Modelica.Icons.Record;
        end ThermodynamicState;

        replaceable partial model BaseProperties  "Base properties (p, d, T, h, u, R, MM and, if applicable, X and Xi) of a medium" 
          InputAbsolutePressure p "Absolute pressure of medium";
          InputMassFraction[nXi] Xi(start = reference_X[1:nXi]) "Structurally independent mass fractions";
          InputSpecificEnthalpy h "Specific enthalpy of medium";
          Density d "Density of medium";
          Temperature T "Temperature of medium";
          MassFraction[nX] X(start = reference_X) "Mass fractions (= (component mass)/total mass  m_i/m)";
          SpecificInternalEnergy u "Specific internal energy of medium";
          SpecificHeatCapacity R "Gas constant (of mixture if applicable)";
          MolarMass MM "Molar mass (of mixture or single fluid)";
          ThermodynamicState state "Thermodynamic state record for optional functions";
          parameter Boolean preferredMediumStates = false "= true if StateSelect.prefer shall be used for the independent property variables of the medium" annotation(Evaluate = true);
          parameter Boolean standardOrderComponents = true "If true, and reducedX = true, the last element of X will be computed from the other ones";
          .Modelica.SIunits.Conversions.NonSIunits.Temperature_degC T_degC = Modelica.SIunits.Conversions.to_degC(T) "Temperature of medium in [degC]";
          .Modelica.SIunits.Conversions.NonSIunits.Pressure_bar p_bar = Modelica.SIunits.Conversions.to_bar(p) "Absolute pressure of medium in [bar]";
          connector InputAbsolutePressure = input .Modelica.SIunits.AbsolutePressure "Pressure as input signal connector";
          connector InputSpecificEnthalpy = input .Modelica.SIunits.SpecificEnthalpy "Specific enthalpy as input signal connector";
          connector InputMassFraction = input .Modelica.SIunits.MassFraction "Mass fraction as input signal connector";
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
              assert(X[i] >= (-1.e-5) and X[i] <= 1 + 1.e-5, "Mass fraction X[" + String(i) + "] = " + String(X[i]) + "of substance " + substanceNames[i] + "\nof medium " + mediumName + " is not in the range 0..1");
            end for;
          end if;
          assert(p >= 0.0, "Pressure (= " + String(p) + " Pa) of medium \"" + mediumName + "\" is negative\n(Temperature = " + String(T) + " K)");
        end BaseProperties;

        replaceable partial function setState_pTX  "Return thermodynamic state as function of p, T and composition X or Xi" 
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input Temperature T "Temperature";
          input MassFraction[:] X = reference_X "Mass fractions";
          output ThermodynamicState state "Thermodynamic state record";
        end setState_pTX;

        replaceable partial function setState_psX  "Return thermodynamic state as function of p, s and composition X or Xi" 
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input SpecificEntropy s "Specific entropy";
          input MassFraction[:] X = reference_X "Mass fractions";
          output ThermodynamicState state "Thermodynamic state record";
        end setState_psX;

        replaceable partial function setSmoothState  "Return thermodynamic state so that it smoothly approximates: if x > 0 then state_a else state_b" 
          extends Modelica.Icons.Function;
          input Real x "m_flow or dp";
          input ThermodynamicState state_a "Thermodynamic state if x > 0";
          input ThermodynamicState state_b "Thermodynamic state if x < 0";
          input Real x_small(min = 0) "Smooth transition in the region -x_small < x < x_small";
          output ThermodynamicState state "Smooth thermodynamic state for all x (continuous and differentiable)";
        end setSmoothState;

        replaceable partial function dynamicViscosity  "Return dynamic viscosity" 
          extends Modelica.Icons.Function;
          input ThermodynamicState state "Thermodynamic state record";
          output DynamicViscosity eta "Dynamic viscosity";
        end dynamicViscosity;

        replaceable partial function thermalConductivity  "Return thermal conductivity" 
          extends Modelica.Icons.Function;
          input ThermodynamicState state "Thermodynamic state record";
          output ThermalConductivity lambda "Thermal conductivity";
        end thermalConductivity;

        replaceable partial function pressure  "Return pressure" 
          extends Modelica.Icons.Function;
          input ThermodynamicState state "Thermodynamic state record";
          output AbsolutePressure p "Pressure";
        end pressure;

        replaceable partial function temperature  "Return temperature" 
          extends Modelica.Icons.Function;
          input ThermodynamicState state "Thermodynamic state record";
          output Temperature T "Temperature";
        end temperature;

        replaceable partial function density  "Return density" 
          extends Modelica.Icons.Function;
          input ThermodynamicState state "Thermodynamic state record";
          output Density d "Density";
        end density;

        replaceable partial function specificEnthalpy  "Return specific enthalpy" 
          extends Modelica.Icons.Function;
          input ThermodynamicState state "Thermodynamic state record";
          output SpecificEnthalpy h "Specific enthalpy";
        end specificEnthalpy;

        replaceable partial function specificInternalEnergy  "Return specific internal energy" 
          extends Modelica.Icons.Function;
          input ThermodynamicState state "Thermodynamic state record";
          output SpecificEnergy u "Specific internal energy";
        end specificInternalEnergy;

        replaceable partial function specificEntropy  "Return specific entropy" 
          extends Modelica.Icons.Function;
          input ThermodynamicState state "Thermodynamic state record";
          output SpecificEntropy s "Specific entropy";
        end specificEntropy;

        replaceable partial function specificGibbsEnergy  "Return specific Gibbs energy" 
          extends Modelica.Icons.Function;
          input ThermodynamicState state "Thermodynamic state record";
          output SpecificEnergy g "Specific Gibbs energy";
        end specificGibbsEnergy;

        replaceable partial function specificHelmholtzEnergy  "Return specific Helmholtz energy" 
          extends Modelica.Icons.Function;
          input ThermodynamicState state "Thermodynamic state record";
          output SpecificEnergy f "Specific Helmholtz energy";
        end specificHelmholtzEnergy;

        replaceable partial function specificHeatCapacityCp  "Return specific heat capacity at constant pressure" 
          extends Modelica.Icons.Function;
          input ThermodynamicState state "Thermodynamic state record";
          output SpecificHeatCapacity cp "Specific heat capacity at constant pressure";
        end specificHeatCapacityCp;

        replaceable partial function specificHeatCapacityCv  "Return specific heat capacity at constant volume" 
          extends Modelica.Icons.Function;
          input ThermodynamicState state "Thermodynamic state record";
          output SpecificHeatCapacity cv "Specific heat capacity at constant volume";
        end specificHeatCapacityCv;

        replaceable partial function isentropicExponent  "Return isentropic exponent" 
          extends Modelica.Icons.Function;
          input ThermodynamicState state "Thermodynamic state record";
          output IsentropicExponent gamma "Isentropic exponent";
        end isentropicExponent;

        replaceable partial function isentropicEnthalpy  "Return isentropic enthalpy" 
          extends Modelica.Icons.Function;
          input AbsolutePressure p_downstream "Downstream pressure";
          input ThermodynamicState refState "Reference state for entropy";
          output SpecificEnthalpy h_is "Isentropic enthalpy";
        end isentropicEnthalpy;

        replaceable partial function velocityOfSound  "Return velocity of sound" 
          extends Modelica.Icons.Function;
          input ThermodynamicState state "Thermodynamic state record";
          output VelocityOfSound a "Velocity of sound";
        end velocityOfSound;

        replaceable partial function isobaricExpansionCoefficient  "Return overall the isobaric expansion coefficient beta" 
          extends Modelica.Icons.Function;
          input ThermodynamicState state "Thermodynamic state record";
          output IsobaricExpansionCoefficient beta "Isobaric expansion coefficient";
        end isobaricExpansionCoefficient;

        replaceable partial function isothermalCompressibility  "Return overall the isothermal compressibility factor" 
          extends Modelica.Icons.Function;
          input ThermodynamicState state "Thermodynamic state record";
          output .Modelica.SIunits.IsothermalCompressibility kappa "Isothermal compressibility";
        end isothermalCompressibility;

        replaceable partial function density_derp_T  "Return density derivative w.r.t. pressure at const temperature" 
          extends Modelica.Icons.Function;
          input ThermodynamicState state "Thermodynamic state record";
          output DerDensityByPressure ddpT "Density derivative w.r.t. pressure";
        end density_derp_T;

        replaceable partial function density_derT_p  "Return density derivative w.r.t. temperature at constant pressure" 
          extends Modelica.Icons.Function;
          input ThermodynamicState state "Thermodynamic state record";
          output DerDensityByTemperature ddTp "Density derivative w.r.t. temperature";
        end density_derT_p;

        replaceable partial function molarMass  "Return the molar mass of the medium" 
          extends Modelica.Icons.Function;
          input ThermodynamicState state "Thermodynamic state record";
          output MolarMass MM "Mixture molar mass";
        end molarMass;

        replaceable function specificEnthalpy_pTX  "Return specific enthalpy from p, T, and X or Xi" 
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input Temperature T "Temperature";
          input MassFraction[:] X = reference_X "Mass fractions";
          output SpecificEnthalpy h "Specific enthalpy";
        algorithm
          h := specificEnthalpy(setState_pTX(p, T, X));
          annotation(inverse(T = temperature_phX(p, h, X))); 
        end specificEnthalpy_pTX;

        replaceable function specificEnthalpy_psX  "Return specific enthalpy from p, s, and X or Xi" 
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input SpecificEntropy s "Specific entropy";
          input MassFraction[:] X = reference_X "Mass fractions";
          output SpecificEnthalpy h "Specific enthalpy";
        algorithm
          h := specificEnthalpy(setState_psX(p, s, X));
        end specificEnthalpy_psX;

        type MassFlowRate = .Modelica.SIunits.MassFlowRate(quantity = "MassFlowRate." + mediumName, min = -1.0e5, max = 1.e5) "Type for mass flow rate with medium specific attributes";
      end PartialMedium;

      partial package PartialMixtureMedium  "Base class for pure substances of several chemical substances" 
        extends PartialMedium(redeclare replaceable record FluidConstants = Modelica.Media.Interfaces.Types.IdealGas.FluidConstants);

        redeclare replaceable record extends ThermodynamicState  "Thermodynamic state variables" 
          AbsolutePressure p "Absolute pressure of medium";
          Temperature T "Temperature of medium";
          MassFraction[nX] X(start = reference_X) "Mass fractions (= (component mass)/total mass  m_i/m)";
        end ThermodynamicState;

        constant FluidConstants[nS] fluidConstants "Constant data for the fluid";

        replaceable function gasConstant  "Return the gas constant of the mixture (also for liquids)" 
          extends Modelica.Icons.Function;
          input ThermodynamicState state "Thermodynamic state";
          output .Modelica.SIunits.SpecificHeatCapacity R "Mixture gas constant";
        end gasConstant;

        function massToMoleFractions  "Return mole fractions from mass fractions X" 
          extends Modelica.Icons.Function;
          input .Modelica.SIunits.MassFraction[:] X "Mass fractions of mixture";
          input .Modelica.SIunits.MolarMass[:] MMX "Molar masses of components";
          output .Modelica.SIunits.MoleFraction[size(X, 1)] moleFractions "Mole fractions of gas mixture";
        protected
          Real[size(X, 1)] invMMX "Inverses of molar weights";
          .Modelica.SIunits.MolarMass Mmix "Molar mass of mixture";
        algorithm
          for i in 1:size(X, 1) loop
            invMMX[i] := 1 / MMX[i];
          end for;
          Mmix := 1 / (X * invMMX);
          for i in 1:size(X, 1) loop
            moleFractions[i] := Mmix * X[i] / MMX[i];
          end for;
          annotation(smoothOrder = 5); 
        end massToMoleFractions;
      end PartialMixtureMedium;

      package Choices  "Types, constants to define menu choices" 
        extends Modelica.Icons.Package;
        type IndependentVariables = enumeration(T "Temperature", pT "Pressure, Temperature", ph "Pressure, Specific Enthalpy", phX "Pressure, Specific Enthalpy, Mass Fraction", pTX "Pressure, Temperature, Mass Fractions", dTX "Density, Temperature, Mass Fractions") "Enumeration defining the independent variables of a medium";
        type ReferenceEnthalpy = enumeration(ZeroAt0K "The enthalpy is 0 at 0 K (default), if the enthalpy of formation is excluded", ZeroAt25C "The enthalpy is 0 at 25 degC, if the enthalpy of formation is excluded", UserDefined "The user-defined reference enthalpy is used at 293.15 K (25 degC)") "Enumeration defining the reference enthalpy of a medium" annotation(Evaluate = true);
      end Choices;

      package Types  "Types to be used in fluid models" 
        extends Modelica.Icons.Package;
        type AbsolutePressure = .Modelica.SIunits.AbsolutePressure(min = 0, max = 1.e8, nominal = 1.e5, start = 1.e5) "Type for absolute pressure with medium specific attributes";
        type Density = .Modelica.SIunits.Density(min = 0, max = 1.e5, nominal = 1, start = 1) "Type for density with medium specific attributes";
        type DynamicViscosity = .Modelica.SIunits.DynamicViscosity(min = 0, max = 1.e8, nominal = 1.e-3, start = 1.e-3) "Type for dynamic viscosity with medium specific attributes";
        type MassFraction = Real(quantity = "MassFraction", final unit = "kg/kg", min = 0, max = 1, nominal = 0.1) "Type for mass fraction with medium specific attributes";
        type MoleFraction = Real(quantity = "MoleFraction", final unit = "mol/mol", min = 0, max = 1, nominal = 0.1) "Type for mole fraction with medium specific attributes";
        type MolarMass = .Modelica.SIunits.MolarMass(min = 0.001, max = 0.25, nominal = 0.032) "Type for molar mass with medium specific attributes";
        type MolarVolume = .Modelica.SIunits.MolarVolume(min = 1e-6, max = 1.0e6, nominal = 1.0) "Type for molar volume with medium specific attributes";
        type IsentropicExponent = .Modelica.SIunits.RatioOfSpecificHeatCapacities(min = 1, max = 500000, nominal = 1.2, start = 1.2) "Type for isentropic exponent with medium specific attributes";
        type SpecificEnergy = .Modelica.SIunits.SpecificEnergy(min = -1.0e8, max = 1.e8, nominal = 1.e6) "Type for specific energy with medium specific attributes";
        type SpecificInternalEnergy = SpecificEnergy "Type for specific internal energy with medium specific attributes";
        type SpecificEnthalpy = .Modelica.SIunits.SpecificEnthalpy(min = -1.0e10, max = 1.e10, nominal = 1.e6) "Type for specific enthalpy with medium specific attributes";
        type SpecificEntropy = .Modelica.SIunits.SpecificEntropy(min = -1.e7, max = 1.e7, nominal = 1.e3) "Type for specific entropy with medium specific attributes";
        type SpecificHeatCapacity = .Modelica.SIunits.SpecificHeatCapacity(min = 0, max = 1.e7, nominal = 1.e3, start = 1.e3) "Type for specific heat capacity with medium specific attributes";
        type Temperature = .Modelica.SIunits.Temperature(min = 1, max = 1.e4, nominal = 300, start = 288.15) "Type for temperature with medium specific attributes";
        type ThermalConductivity = .Modelica.SIunits.ThermalConductivity(min = 0, max = 500, nominal = 1, start = 1) "Type for thermal conductivity with medium specific attributes";
        type VelocityOfSound = .Modelica.SIunits.Velocity(min = 0, max = 1.e5, nominal = 1000, start = 1000) "Type for velocity of sound with medium specific attributes";
        type ExtraProperty = Real(min = 0.0, start = 1.0) "Type for unspecified, mass-specific property transported by flow";
        type IsobaricExpansionCoefficient = Real(min = 0, max = 1.0e8, unit = "1/K") "Type for isobaric expansion coefficient with medium specific attributes";
        type DipoleMoment = Real(min = 0.0, max = 2.0, unit = "debye", quantity = "ElectricDipoleMoment") "Type for dipole moment with medium specific attributes";
        type DerDensityByPressure = .Modelica.SIunits.DerDensityByPressure "Type for partial derivative of density with respect to pressure with medium specific attributes";
        type DerDensityByTemperature = .Modelica.SIunits.DerDensityByTemperature "Type for partial derivative of density with respect to temperature with medium specific attributes";

        package Basic  "The most basic version of a record used in several degrees of detail" 
          extends Icons.Package;

          record FluidConstants  "Critical, triple, molecular and other standard data of fluid" 
            extends Modelica.Icons.Record;
            String iupacName "Complete IUPAC name (or common name, if non-existent)";
            String casRegistryNumber "Chemical abstracts sequencing number (if it exists)";
            String chemicalFormula "Chemical formula, (brutto, nomenclature according to Hill";
            String structureFormula "Chemical structure formula";
            MolarMass molarMass "Molar mass";
          end FluidConstants;
        end Basic;

        package IdealGas  "The ideal gas version of a record used in several degrees of detail" 
          extends Icons.Package;

          record FluidConstants  "Extended fluid constants" 
            extends Modelica.Media.Interfaces.Types.Basic.FluidConstants;
            Temperature criticalTemperature "Critical temperature";
            AbsolutePressure criticalPressure "Critical pressure";
            MolarVolume criticalMolarVolume "Critical molar Volume";
            Real acentricFactor "Pitzer acentric factor";
            Temperature meltingPoint "Melting point at 101325 Pa";
            Temperature normalBoilingPoint "Normal boiling point (at 101325 Pa)";
            DipoleMoment dipoleMoment "Dipole moment of molecule in Debye (1 debye = 3.33564e10-30 C.m)";
            Boolean hasIdealGasHeatCapacity = false "True if ideal gas heat capacity is available";
            Boolean hasCriticalData = false "True if critical data are known";
            Boolean hasDipoleMoment = false "True if a dipole moment known";
            Boolean hasFundamentalEquation = false "True if a fundamental equation";
            Boolean hasLiquidHeatCapacity = false "True if liquid heat capacity is available";
            Boolean hasSolidHeatCapacity = false "True if solid heat capacity is available";
            Boolean hasAccurateViscosityData = false "True if accurate data for a viscosity function is available";
            Boolean hasAccurateConductivityData = false "True if accurate data for thermal conductivity is available";
            Boolean hasVapourPressureCurve = false "True if vapour pressure data, e.g., Antoine coefficients are known";
            Boolean hasAcentricFactor = false "True if Pitzer accentric factor is known";
            SpecificEnthalpy HCRIT0 = 0.0 "Critical specific enthalpy of the fundamental equation";
            SpecificEntropy SCRIT0 = 0.0 "Critical specific entropy of the fundamental equation";
            SpecificEnthalpy deltah = 0.0 "Difference between specific enthalpy model (h_m) and f.eq. (h_f) (h_m - h_f)";
            SpecificEntropy deltas = 0.0 "Difference between specific enthalpy model (s_m) and f.eq. (s_f) (s_m - s_f)";
          end FluidConstants;
        end IdealGas;
      end Types;
    end Interfaces;

    package Common  "Data structures and fundamental functions for fluid properties" 
      extends Modelica.Icons.Package;
      constant Real MINPOS = 1.0e-9 "Minimal value for physical variables which are always > 0.0";

      function smoothStep  "Approximation of a general step, such that the characteristic is continuous and differentiable" 
        extends Modelica.Icons.Function;
        input Real x "Abscissa value";
        input Real y1 "Ordinate value for x > 0";
        input Real y2 "Ordinate value for x < 0";
        input Real x_small(min = 0) = 1e-5 "Approximation of step for -x_small <= x <= x_small; x_small > 0 required";
        output Real y "Ordinate value to approximate y = if x > 0 then y1 else y2";
      algorithm
        y := smooth(1, if x > x_small then y1 else if x < (-x_small) then y2 else if abs(x_small) > 0 then x / x_small * ((x / x_small) ^ 2 - 3) * (y2 - y1) / 4 + (y1 + y2) / 2 else (y1 + y2) / 2);
        annotation(Inline = true, smoothOrder = 1); 
      end smoothStep;

      package OneNonLinearEquation  "Determine solution of a non-linear algebraic equation in one unknown without derivatives in a reliable and efficient way" 
        extends Modelica.Icons.Package;

        replaceable record f_nonlinear_Data  "Data specific for function f_nonlinear" 
          extends Modelica.Icons.Record;
        end f_nonlinear_Data;

        replaceable partial function f_nonlinear  "Nonlinear algebraic equation in one unknown: y = f_nonlinear(x,p,X)" 
          extends Modelica.Icons.Function;
          input Real x "Independent variable of function";
          input Real p = 0.0 "Disregarded variables (here always used for pressure)";
          input Real[:] X = fill(0, 0) "Disregarded variables (her always used for composition)";
          input f_nonlinear_Data f_nonlinear_data "Additional data for the function";
          output Real y "= f_nonlinear(x)";
        end f_nonlinear;

        replaceable function solve  "Solve f_nonlinear(x_zero)=y_zero; f_nonlinear(x_min) - y_zero and f_nonlinear(x_max)-y_zero must have different sign" 
          extends Modelica.Icons.Function;
          input Real y_zero "Determine x_zero, such that f_nonlinear(x_zero) = y_zero";
          input Real x_min "Minimum value of x";
          input Real x_max "Maximum value of x";
          input Real pressure = 0.0 "Disregarded variables (here always used for pressure)";
          input Real[:] X = fill(0, 0) "Disregarded variables (here always used for composition)";
          input f_nonlinear_Data f_nonlinear_data "Additional data for function f_nonlinear";
          input Real x_tol = 100 * Modelica.Constants.eps "Relative tolerance of the result";
          output Real x_zero "f_nonlinear(x_zero) = y_zero";
        protected
          constant Real eps = Modelica.Constants.eps "Machine epsilon";
          constant Real x_eps = 1e-10 "Slight modification of x_min, x_max, since x_min, x_max are usually exactly at the borders T_min/h_min and then small numeric noise may make the interval invalid";
          Real x_min2 = x_min - x_eps;
          Real x_max2 = x_max + x_eps;
          Real a = x_min2 "Current best minimum interval value";
          Real b = x_max2 "Current best maximum interval value";
          Real c "Intermediate point a <= c <= b";
          Real d;
          Real e "b - a";
          Real m;
          Real s;
          Real p;
          Real q;
          Real r;
          Real tol;
          Real fa "= f_nonlinear(a) - y_zero";
          Real fb "= f_nonlinear(b) - y_zero";
          Real fc;
          Boolean found = false;
        algorithm
          fa := f_nonlinear(x_min2, pressure, X, f_nonlinear_data) - y_zero;
          fb := f_nonlinear(x_max2, pressure, X, f_nonlinear_data) - y_zero;
          fc := fb;
          if fa > 0.0 and fb > 0.0 or fa < 0.0 and fb < 0.0 then
            .Modelica.Utilities.Streams.error("The arguments x_min and x_max to OneNonLinearEquation.solve(..)\n" + "do not bracket the root of the single non-linear equation:\n" + "  x_min  = " + String(x_min2) + "\n" + "  x_max  = " + String(x_max2) + "\n" + "  y_zero = " + String(y_zero) + "\n" + "  fa = f(x_min) - y_zero = " + String(fa) + "\n" + "  fb = f(x_max) - y_zero = " + String(fb) + "\n" + "fa and fb must have opposite sign which is not the case");
          else
          end if;
          c := a;
          fc := fa;
          e := b - a;
          d := e;
          while not found loop
            if abs(fc) < abs(fb) then
              a := b;
              b := c;
              c := a;
              fa := fb;
              fb := fc;
              fc := fa;
            else
            end if;
            tol := 2 * eps * abs(b) + x_tol;
            m := (c - b) / 2;
            if abs(m) <= tol or fb == 0.0 then
              found := true;
              x_zero := b;
            else
              if abs(e) < tol or abs(fa) <= abs(fb) then
                e := m;
                d := e;
              else
                s := fb / fa;
                if a == c then
                  p := 2 * m * s;
                  q := 1 - s;
                else
                  q := fa / fc;
                  r := fb / fc;
                  p := s * (2 * m * q * (q - r) - (b - a) * (r - 1));
                  q := (q - 1) * (r - 1) * (s - 1);
                end if;
                if p > 0 then
                  q := -q;
                else
                  p := -p;
                end if;
                s := e;
                e := d;
                if 2 * p < 3 * m * q - abs(tol * q) and p < abs(0.5 * s * q) then
                  d := p / q;
                else
                  e := m;
                  d := e;
                end if;
              end if;
              a := b;
              fa := fb;
              b := b + (if abs(d) > tol then d else if m > 0 then tol else -tol);
              fb := f_nonlinear(b, pressure, X, f_nonlinear_data) - y_zero;
              if fb > 0 and fc > 0 or fb < 0 and fc < 0 then
                c := a;
                fc := fa;
                e := b - a;
                d := e;
              else
              end if;
            end if;
          end while;
        end solve;
      end OneNonLinearEquation;
    end Common;

    package IdealGases  "Data and models of ideal gases (single, fixed and dynamic mixtures) from NASA source" 
      extends Modelica.Icons.VariantsPackage;

      package Common  "Common packages and data for the ideal gas models" 
        extends Modelica.Icons.Package;

        record DataRecord  "Coefficient data record for properties of ideal gases based on NASA source" 
          extends Modelica.Icons.Record;
          String name "Name of ideal gas";
          .Modelica.SIunits.MolarMass MM "Molar mass";
          .Modelica.SIunits.SpecificEnthalpy Hf "Enthalpy of formation at 298.15K";
          .Modelica.SIunits.SpecificEnthalpy H0 "H0(298.15K) - H0(0K)";
          .Modelica.SIunits.Temperature Tlimit "Temperature limit between low and high data sets";
          Real[7] alow "Low temperature coefficients a";
          Real[2] blow "Low temperature constants b";
          Real[7] ahigh "High temperature coefficients a";
          Real[2] bhigh "High temperature constants b";
          .Modelica.SIunits.SpecificHeatCapacity R "Gas constant";
        end DataRecord;

        package Functions  "Basic Functions for ideal gases: cp, h, s, thermal conductivity, viscosity" 
          extends Modelica.Icons.Package;
          constant Boolean excludeEnthalpyOfFormation = true "If true, enthalpy of formation Hf is not included in specific enthalpy h";
          constant Modelica.Media.Interfaces.Choices.ReferenceEnthalpy referenceChoice = Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt0K "Choice of reference enthalpy";
          constant Modelica.Media.Interfaces.Types.SpecificEnthalpy h_offset = 0.0 "User defined offset for reference enthalpy, if referenceChoice = UserDefined";

          function cp_T  "Compute specific heat capacity at constant pressure from temperature and gas data" 
            extends Modelica.Icons.Function;
            input IdealGases.Common.DataRecord data "Ideal gas data";
            input .Modelica.SIunits.Temperature T "Temperature";
            output .Modelica.SIunits.SpecificHeatCapacity cp "Specific heat capacity at temperature T";
          algorithm
            cp := smooth(0, if T < data.Tlimit then data.R * (1 / (T * T) * (data.alow[1] + T * (data.alow[2] + T * (1. * data.alow[3] + T * (data.alow[4] + T * (data.alow[5] + T * (data.alow[6] + data.alow[7] * T))))))) else data.R * (1 / (T * T) * (data.ahigh[1] + T * (data.ahigh[2] + T * (1. * data.ahigh[3] + T * (data.ahigh[4] + T * (data.ahigh[5] + T * (data.ahigh[6] + data.ahigh[7] * T))))))));
            annotation(Inline = true, smoothOrder = 2); 
          end cp_T;

          function h_T  "Compute specific enthalpy from temperature and gas data; reference is decided by the
              refChoice input, or by the referenceChoice package constant by default" 
            extends Modelica.Icons.Function;
            input IdealGases.Common.DataRecord data "Ideal gas data";
            input .Modelica.SIunits.Temperature T "Temperature";
            input Boolean exclEnthForm = excludeEnthalpyOfFormation "If true, enthalpy of formation Hf is not included in specific enthalpy h";
            input Modelica.Media.Interfaces.Choices.ReferenceEnthalpy refChoice = referenceChoice "Choice of reference enthalpy";
            input .Modelica.SIunits.SpecificEnthalpy h_off = h_offset "User defined offset for reference enthalpy, if referenceChoice = UserDefined";
            output .Modelica.SIunits.SpecificEnthalpy h "Specific enthalpy at temperature T";
          algorithm
            h := smooth(0, (if T < data.Tlimit then data.R * (((-data.alow[1]) + T * (data.blow[1] + data.alow[2] * Math.log(T) + T * (1. * data.alow[3] + T * (0.5 * data.alow[4] + T * (1 / 3 * data.alow[5] + T * (0.25 * data.alow[6] + 0.2 * data.alow[7] * T)))))) / T) else data.R * (((-data.ahigh[1]) + T * (data.bhigh[1] + data.ahigh[2] * Math.log(T) + T * (1. * data.ahigh[3] + T * (0.5 * data.ahigh[4] + T * (1 / 3 * data.ahigh[5] + T * (0.25 * data.ahigh[6] + 0.2 * data.ahigh[7] * T)))))) / T)) + (if exclEnthForm then -data.Hf else 0.0) + (if refChoice == .Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt0K then data.H0 else 0.0) + (if refChoice == .Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.UserDefined then h_off else 0.0));
            annotation(Inline = false, smoothOrder = 2); 
          end h_T;

          function s0_T  "Compute specific entropy from temperature and gas data" 
            extends Modelica.Icons.Function;
            input IdealGases.Common.DataRecord data "Ideal gas data";
            input .Modelica.SIunits.Temperature T "Temperature";
            output .Modelica.SIunits.SpecificEntropy s "Specific entropy at temperature T";
          algorithm
            s := if T < data.Tlimit then data.R * (data.blow[2] - 0.5 * data.alow[1] / (T * T) - data.alow[2] / T + data.alow[3] * Math.log(T) + T * (data.alow[4] + T * (0.5 * data.alow[5] + T * (1 / 3 * data.alow[6] + 0.25 * data.alow[7] * T)))) else data.R * (data.bhigh[2] - 0.5 * data.ahigh[1] / (T * T) - data.ahigh[2] / T + data.ahigh[3] * Math.log(T) + T * (data.ahigh[4] + T * (0.5 * data.ahigh[5] + T * (1 / 3 * data.ahigh[6] + 0.25 * data.ahigh[7] * T))));
            annotation(Inline = true, smoothOrder = 2); 
          end s0_T;

          function dynamicViscosityLowPressure  "Dynamic viscosity of low pressure gases" 
            extends Modelica.Icons.Function;
            input .Modelica.SIunits.Temp_K T "Gas temperature";
            input .Modelica.SIunits.Temp_K Tc "Critical temperature of gas";
            input .Modelica.SIunits.MolarMass M "Molar mass of gas";
            input .Modelica.SIunits.MolarVolume Vc "Critical molar volume of gas";
            input Real w "Acentric factor of gas";
            input Modelica.Media.Interfaces.Types.DipoleMoment mu "Dipole moment of gas molecule";
            input Real k = 0.0 "Special correction for highly polar substances";
            output .Modelica.SIunits.DynamicViscosity eta "Dynamic viscosity of gas";
          protected
            parameter Real Const1_SI = 40.785 * 10 ^ (-9.5) "Constant in formula for eta converted to SI units";
            parameter Real Const2_SI = 131.3 / 1000.0 "Constant in formula for mur converted to SI units";
            Real mur = Const2_SI * mu / sqrt(Vc * Tc) "Dimensionless dipole moment of gas molecule";
            Real Fc = 1 - 0.2756 * w + 0.059035 * mur ^ 4 + k "Factor to account for molecular shape and polarities of gas";
            Real Tstar "Dimensionless temperature defined by equation below";
            Real Ov "Viscosity collision integral for the gas";
          algorithm
            Tstar := 1.2593 * T / Tc;
            Ov := 1.16145 * Tstar ^ (-0.14874) + 0.52487 * Modelica.Math.exp(-0.7732 * Tstar) + 2.16178 * Modelica.Math.exp(-2.43787 * Tstar);
            eta := Const1_SI * Fc * sqrt(M * T) / (Vc ^ (2 / 3) * Ov);
            annotation(smoothOrder = 2); 
          end dynamicViscosityLowPressure;

          function thermalConductivityEstimate  "Thermal conductivity of polyatomic gases(Eucken and Modified Eucken correlation)" 
            extends Modelica.Icons.Function;
            input Modelica.Media.Interfaces.Types.SpecificHeatCapacity Cp "Constant pressure heat capacity";
            input Modelica.Media.Interfaces.Types.DynamicViscosity eta "Dynamic viscosity";
            input Integer method(min = 1, max = 2) = 1 "1: Eucken Method, 2: Modified Eucken Method";
            input IdealGases.Common.DataRecord data "Ideal gas data";
            output Modelica.Media.Interfaces.Types.ThermalConductivity lambda "Thermal conductivity [W/(m.k)]";
          algorithm
            lambda := if method == 1 then eta * (Cp - data.R + 9 / 4 * data.R) else eta * (Cp - data.R) * (1.32 + 1.77 / (Cp / Modelica.Constants.R - 1.0));
            annotation(smoothOrder = 2); 
          end thermalConductivityEstimate;
        end Functions;

        partial package MixtureGasNasa  "Medium model of a mixture of ideal gases based on NASA source" 
          extends Modelica.Media.Interfaces.PartialMixtureMedium(ThermoStates = Modelica.Media.Interfaces.Choices.IndependentVariables.pTX, substanceNames = data[:].name, reducedX = false, singleState = false, reference_X = fill(1 / nX, nX), SpecificEnthalpy(start = if referenceChoice == .Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt0K then 3e5 else if referenceChoice == .Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.UserDefined then h_offset else 0, nominal = 1.0e5), Density(start = 10, nominal = 10), AbsolutePressure(start = 10e5, nominal = 10e5), Temperature(min = 200, max = 6000, start = 500, nominal = 500));

          redeclare record extends ThermodynamicState  "Thermodynamic state variables" end ThermodynamicState;

          constant Modelica.Media.IdealGases.Common.DataRecord[:] data "Data records of ideal gas substances";
          constant Boolean excludeEnthalpyOfFormation = true "If true, enthalpy of formation Hf is not included in specific enthalpy h";
          constant .Modelica.Media.Interfaces.Choices.ReferenceEnthalpy referenceChoice = .Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt0K "Choice of reference enthalpy";
          constant SpecificEnthalpy h_offset = 0.0 "User defined offset for reference enthalpy, if referenceChoice = UserDefined";
          constant MolarMass[nX] MMX = data[:].MM "Molar masses of components";
          constant Integer methodForThermalConductivity(min = 1, max = 2) = 1;

          redeclare replaceable model extends BaseProperties(T(stateSelect = if preferredMediumStates then StateSelect.prefer else StateSelect.default), p(stateSelect = if preferredMediumStates then StateSelect.prefer else StateSelect.default), Xi(each stateSelect = if preferredMediumStates then StateSelect.prefer else StateSelect.default), final standardOrderComponents = true)  "Base properties (p, d, T, h, u, R, MM, X, and Xi of NASA mixture gas" 
          equation
            assert(T >= 200 and T <= 6000, "
          Temperature T (=" + String(T) + " K = 200 K) is not in the allowed range
          200 K <= T <= 6000 K
          required from medium model \"" + mediumName + "\".");
            MM = molarMass(state);
            h = h_TX(T, X);
            R = data.R * X;
            u = h - R * T;
            d = p / (R * T);
            state.T = T;
            state.p = p;
            state.X = if fixedX then reference_X else X;
          end BaseProperties;

          redeclare function setState_pTX  "Return thermodynamic state as function of p, T and composition X" 
            extends Modelica.Icons.Function;
            input AbsolutePressure p "Pressure";
            input Temperature T "Temperature";
            input MassFraction[:] X = reference_X "Mass fractions";
            output ThermodynamicState state;
          algorithm
            state := if size(X, 1) == 0 then ThermodynamicState(p = p, T = T, X = reference_X) else if size(X, 1) == nX then ThermodynamicState(p = p, T = T, X = X) else ThermodynamicState(p = p, T = T, X = cat(1, X, {1 - sum(X)}));
            annotation(Inline = true, smoothOrder = 2); 
          end setState_pTX;

          redeclare function setState_psX  "Return thermodynamic state as function of p, s and composition X" 
            extends Modelica.Icons.Function;
            input AbsolutePressure p "Pressure";
            input SpecificEntropy s "Specific entropy";
            input MassFraction[:] X = reference_X "Mass fractions";
            output ThermodynamicState state;
          algorithm
            state := if size(X, 1) == 0 then ThermodynamicState(p = p, T = T_psX(p, s, reference_X), X = reference_X) else if size(X, 1) == nX then ThermodynamicState(p = p, T = T_psX(p, s, X), X = X) else ThermodynamicState(p = p, T = T_psX(p, s, X), X = cat(1, X, {1 - sum(X)}));
            annotation(Inline = true, smoothOrder = 2); 
          end setState_psX;

          redeclare function extends setSmoothState  "Return thermodynamic state so that it smoothly approximates: if x > 0 then state_a else state_b" 
          algorithm
            state := ThermodynamicState(p = Media.Common.smoothStep(x, state_a.p, state_b.p, x_small), T = Media.Common.smoothStep(x, state_a.T, state_b.T, x_small), X = Media.Common.smoothStep(x, state_a.X, state_b.X, x_small));
            annotation(Inline = true, smoothOrder = 2); 
          end setSmoothState;

          redeclare function extends pressure  "Return pressure of ideal gas" 
          algorithm
            p := state.p;
            annotation(Inline = true, smoothOrder = 2); 
          end pressure;

          redeclare function extends temperature  "Return temperature of ideal gas" 
          algorithm
            T := state.T;
            annotation(Inline = true, smoothOrder = 2); 
          end temperature;

          redeclare function extends density  "Return density of ideal gas" 
          algorithm
            d := state.p / (state.X * data.R * state.T);
            annotation(Inline = true, smoothOrder = 3); 
          end density;

          redeclare function extends specificEnthalpy  "Return specific enthalpy" 
            extends Modelica.Icons.Function;
          algorithm
            h := h_TX(state.T, state.X);
            annotation(Inline = true, smoothOrder = 2); 
          end specificEnthalpy;

          redeclare function extends specificInternalEnergy  "Return specific internal energy" 
            extends Modelica.Icons.Function;
          algorithm
            u := h_TX(state.T, state.X) - gasConstant(state) * state.T;
            annotation(Inline = true, smoothOrder = 2); 
          end specificInternalEnergy;

          redeclare function extends specificEntropy  "Return specific entropy" 
          protected
            Real[nX] Y(unit = "mol/mol") = massToMoleFractions(state.X, data.MM) "Molar fractions";
          algorithm
            s := s_TX(state.T, state.X) - sum(state.X[i] * Modelica.Constants.R / MMX[i] * (if state.X[i] < Modelica.Constants.eps then Y[i] else Modelica.Math.log(Y[i] * state.p / reference_p)) for i in 1:nX);
            annotation(Inline = true, smoothOrder = 2); 
          end specificEntropy;

          redeclare function extends specificGibbsEnergy  "Return specific Gibbs energy" 
            extends Modelica.Icons.Function;
          algorithm
            g := h_TX(state.T, state.X) - state.T * specificEntropy(state);
            annotation(Inline = true, smoothOrder = 2); 
          end specificGibbsEnergy;

          redeclare function extends specificHelmholtzEnergy  "Return specific Helmholtz energy" 
            extends Modelica.Icons.Function;
          algorithm
            f := h_TX(state.T, state.X) - gasConstant(state) * state.T - state.T * specificEntropy(state);
            annotation(Inline = true, smoothOrder = 2); 
          end specificHelmholtzEnergy;

          function h_TX  "Return specific enthalpy" 
            extends Modelica.Icons.Function;
            input .Modelica.SIunits.Temperature T "Temperature";
            input MassFraction[nX] X = reference_X "Independent Mass fractions of gas mixture";
            input Boolean exclEnthForm = excludeEnthalpyOfFormation "If true, enthalpy of formation Hf is not included in specific enthalpy h";
            input Modelica.Media.Interfaces.Choices.ReferenceEnthalpy refChoice = referenceChoice "Choice of reference enthalpy";
            input .Modelica.SIunits.SpecificEnthalpy h_off = h_offset "User defined offset for reference enthalpy, if referenceChoice = UserDefined";
            output .Modelica.SIunits.SpecificEnthalpy h "Specific enthalpy at temperature T";
          algorithm
            h := (if fixedX then reference_X else X) * {Modelica.Media.IdealGases.Common.Functions.h_T(data[i], T, exclEnthForm, refChoice, h_off) for i in 1:nX};
            annotation(Inline = false, smoothOrder = 2); 
          end h_TX;

          redeclare function extends gasConstant  "Return gasConstant" 
          algorithm
            R := data.R * state.X;
            annotation(Inline = true, smoothOrder = 3); 
          end gasConstant;

          redeclare function extends specificHeatCapacityCp  "Return specific heat capacity at constant pressure" 
          algorithm
            cp := {Modelica.Media.IdealGases.Common.Functions.cp_T(data[i], state.T) for i in 1:nX} * state.X;
            annotation(Inline = true, smoothOrder = 1); 
          end specificHeatCapacityCp;

          redeclare function extends specificHeatCapacityCv  "Return specific heat capacity at constant volume from temperature and gas data" 
          algorithm
            cv := {Modelica.Media.IdealGases.Common.Functions.cp_T(data[i], state.T) for i in 1:nX} * state.X - data.R * state.X;
            annotation(Inline = true, smoothOrder = 1); 
          end specificHeatCapacityCv;

          function s_TX  "Return temperature dependent part of the entropy, expects full entropy vector" 
            extends Modelica.Icons.Function;
            input Temperature T "Temperature";
            input MassFraction[nX] X "Mass fraction";
            output SpecificEntropy s "Specific entropy";
          algorithm
            s := sum(Modelica.Media.IdealGases.Common.Functions.s0_T(data[i], T) * X[i] for i in 1:size(X, 1));
            annotation(Inline = true, smoothOrder = 2); 
          end s_TX;

          redeclare function extends isentropicExponent  "Return isentropic exponent" 
          algorithm
            gamma := specificHeatCapacityCp(state) / specificHeatCapacityCv(state);
            annotation(Inline = true, smoothOrder = 2); 
          end isentropicExponent;

          redeclare function extends velocityOfSound  "Return velocity of sound" 
            extends Modelica.Icons.Function;
            input ThermodynamicState state "Properties at upstream location";
          algorithm
            a := sqrt(max(0, gasConstant(state) * state.T * specificHeatCapacityCp(state) / specificHeatCapacityCv(state)));
            annotation(Inline = true, smoothOrder = 2); 
          end velocityOfSound;

          function isentropicEnthalpyApproximation  "Approximate method of calculating h_is from upstream properties and downstream pressure" 
            extends Modelica.Icons.Function;
            input AbsolutePressure p2 "Downstream pressure";
            input ThermodynamicState state "Thermodynamic state at upstream location";
            output SpecificEnthalpy h_is "Isentropic enthalpy";
          protected
            SpecificEnthalpy h "Specific enthalpy at upstream location";
            SpecificEnthalpy[nX] h_component "Specific enthalpy at upstream location";
            IsentropicExponent gamma = isentropicExponent(state) "Isentropic exponent";
            MassFraction[nX] X "Complete X-vector";
          algorithm
            X := if reducedX then cat(1, state.X, {1 - sum(state.X)}) else state.X;
            h_component := {Modelica.Media.IdealGases.Common.Functions.h_T(data[i], state.T, excludeEnthalpyOfFormation, referenceChoice, h_offset) for i in 1:nX};
            h := h_component * X;
            h_is := h + gamma / (gamma - 1.0) * (state.T * gasConstant(state)) * ((p2 / state.p) ^ ((gamma - 1) / gamma) - 1.0);
            annotation(smoothOrder = 2); 
          end isentropicEnthalpyApproximation;

          redeclare function extends isentropicEnthalpy  "Return isentropic enthalpy" 
            input Boolean exact = false "Flag whether exact or approximate version should be used";
          algorithm
            h_is := if exact then specificEnthalpy_psX(p_downstream, specificEntropy(refState), refState.X) else isentropicEnthalpyApproximation(p_downstream, refState);
            annotation(Inline = true, smoothOrder = 2); 
          end isentropicEnthalpy;

          function gasMixtureViscosity  "Return viscosities of gas mixtures at low pressures (Wilke method)" 
            extends Modelica.Icons.Function;
            input MoleFraction[:] yi "Mole fractions";
            input MolarMass[size(yi, 1)] M "Mole masses";
            input DynamicViscosity[size(yi, 1)] eta "Pure component viscosities";
            output DynamicViscosity etam "Viscosity of the mixture";
          protected
            Real[size(yi, 1), size(yi, 1)] fi;
          algorithm
            for i in 1:size(eta, 1) loop
              assert(fluidConstants[i].hasDipoleMoment, "Dipole moment for " + fluidConstants[i].chemicalFormula + " not known. Can not compute viscosity.");
              assert(fluidConstants[i].hasCriticalData, "Critical data for " + fluidConstants[i].chemicalFormula + " not known. Can not compute viscosity.");
              for j in 1:size(eta, 1) loop
                if i == 1 then
                  fi[i, j] := (1 + (eta[i] / eta[j]) ^ (1 / 2) * (M[j] / M[i]) ^ (1 / 4)) ^ 2 / (8 * (1 + M[i] / M[j])) ^ (1 / 2);
                elseif j < i then
                  fi[i, j] := eta[i] / eta[j] * M[j] / M[i] * fi[j, i];
                else
                  fi[i, j] := (1 + (eta[i] / eta[j]) ^ (1 / 2) * (M[j] / M[i]) ^ (1 / 4)) ^ 2 / (8 * (1 + M[i] / M[j])) ^ (1 / 2);
                end if;
              end for;
            end for;
            etam := sum(yi[i] * eta[i] / sum(yi[j] * fi[i, j] for j in 1:size(eta, 1)) for i in 1:size(eta, 1));
            annotation(smoothOrder = 2); 
          end gasMixtureViscosity;

          redeclare replaceable function extends dynamicViscosity  "Return mixture dynamic viscosity" 
          protected
            DynamicViscosity[nX] etaX "Component dynamic viscosities";
          algorithm
            for i in 1:nX loop
              etaX[i] := Modelica.Media.IdealGases.Common.Functions.dynamicViscosityLowPressure(state.T, fluidConstants[i].criticalTemperature, fluidConstants[i].molarMass, fluidConstants[i].criticalMolarVolume, fluidConstants[i].acentricFactor, fluidConstants[i].dipoleMoment);
            end for;
            eta := gasMixtureViscosity(massToMoleFractions(state.X, fluidConstants[:].molarMass), fluidConstants[:].molarMass, etaX);
            annotation(smoothOrder = 2); 
          end dynamicViscosity;

          function lowPressureThermalConductivity  "Return thermal conductivities of low-pressure gas mixtures (Mason and Saxena Modification)" 
            extends Modelica.Icons.Function;
            input MoleFraction[:] y "Mole fraction of the components in the gas mixture";
            input Temperature T "Temperature";
            input Temperature[size(y, 1)] Tc "Critical temperatures";
            input AbsolutePressure[size(y, 1)] Pc "Critical pressures";
            input MolarMass[size(y, 1)] M "Molecular weights";
            input ThermalConductivity[size(y, 1)] lambda "Thermal conductivities of the pure gases";
            output ThermalConductivity lambdam "Thermal conductivity of the gas mixture";
          protected
            MolarMass[size(y, 1)] gamma;
            Real[size(y, 1)] Tr "Reduced temperature";
            Real[size(y, 1), size(y, 1)] A "Mason and Saxena Modification";
            constant Real epsilon = 1.0 "Numerical constant near unity";
          algorithm
            for i in 1:size(y, 1) loop
              gamma[i] := 210 * (Tc[i] * M[i] ^ 3 / Pc[i] ^ 4) ^ (1 / 6);
              Tr[i] := T / Tc[i];
            end for;
            for i in 1:size(y, 1) loop
              for j in 1:size(y, 1) loop
                A[i, j] := epsilon * (1 + (gamma[j] * (.Modelica.Math.exp(0.0464 * Tr[i]) - .Modelica.Math.exp(-0.2412 * Tr[i])) / (gamma[i] * (.Modelica.Math.exp(0.0464 * Tr[j]) - .Modelica.Math.exp(-0.2412 * Tr[j])))) ^ (1 / 2) * (M[i] / M[j]) ^ (1 / 4)) ^ 2 / (8 * (1 + M[i] / M[j])) ^ (1 / 2);
              end for;
            end for;
            lambdam := sum(y[i] * lambda[i] / sum(y[j] * A[i, j] for j in 1:size(y, 1)) for i in 1:size(y, 1));
            annotation(smoothOrder = 2); 
          end lowPressureThermalConductivity;

          redeclare replaceable function extends thermalConductivity  "Return thermal conductivity for low pressure gas mixtures" 
            input Integer method = methodForThermalConductivity "Method to compute single component thermal conductivity";
          protected
            ThermalConductivity[nX] lambdaX "Component thermal conductivities";
            DynamicViscosity[nX] eta "Component thermal dynamic viscosities";
            SpecificHeatCapacity[nX] cp "Component heat capacity";
          algorithm
            for i in 1:nX loop
              assert(fluidConstants[i].hasCriticalData, "Critical data for " + fluidConstants[i].chemicalFormula + " not known. Can not compute thermal conductivity.");
              eta[i] := Modelica.Media.IdealGases.Common.Functions.dynamicViscosityLowPressure(state.T, fluidConstants[i].criticalTemperature, fluidConstants[i].molarMass, fluidConstants[i].criticalMolarVolume, fluidConstants[i].acentricFactor, fluidConstants[i].dipoleMoment);
              cp[i] := Modelica.Media.IdealGases.Common.Functions.cp_T(data[i], state.T);
              lambdaX[i] := Modelica.Media.IdealGases.Common.Functions.thermalConductivityEstimate(Cp = cp[i], eta = eta[i], method = method, data = data[i]);
            end for;
            lambda := lowPressureThermalConductivity(massToMoleFractions(state.X, fluidConstants[:].molarMass), state.T, fluidConstants[:].criticalTemperature, fluidConstants[:].criticalPressure, fluidConstants[:].molarMass, lambdaX);
            annotation(smoothOrder = 2); 
          end thermalConductivity;

          redeclare function extends isobaricExpansionCoefficient  "Return isobaric expansion coefficient beta" 
          algorithm
            beta := 1 / state.T;
            annotation(Inline = true, smoothOrder = 2); 
          end isobaricExpansionCoefficient;

          redeclare function extends isothermalCompressibility  "Return isothermal compressibility factor" 
          algorithm
            kappa := 1.0 / state.p;
            annotation(Inline = true, smoothOrder = 2); 
          end isothermalCompressibility;

          redeclare function extends density_derp_T  "Return density derivative by pressure at constant temperature" 
          algorithm
            ddpT := 1 / (state.T * gasConstant(state));
            annotation(Inline = true, smoothOrder = 2); 
          end density_derp_T;

          redeclare function extends density_derT_p  "Return density derivative by temperature at constant pressure" 
          algorithm
            ddTp := -state.p / (state.T * state.T * gasConstant(state));
            annotation(Inline = true, smoothOrder = 2); 
          end density_derT_p;

          redeclare function extends molarMass  "Return molar mass of mixture" 
          algorithm
            MM := 1 / sum(state.X[j] / data[j].MM for j in 1:size(state.X, 1));
            annotation(Inline = true, smoothOrder = 2); 
          end molarMass;

          function T_psX  "Return temperature from pressure, specific entropy and mass fraction" 
            extends Modelica.Icons.Function;
            input AbsolutePressure p "Pressure";
            input SpecificEntropy s "Specific entropy";
            input MassFraction[nX] X "Mass fractions of composition";
            output Temperature T "Temperature";
          protected
            MassFraction[nX] Xfull = if size(X, 1) == nX then X else cat(1, X, {1 - sum(X)});

            package Internal  "Solve h(data,T) for T with given h (use only indirectly via temperature_phX)" 
              extends Modelica.Media.Common.OneNonLinearEquation;

              redeclare record extends f_nonlinear_Data  "Data to be passed to non-linear function" 
                extends Modelica.Media.IdealGases.Common.DataRecord;
              end f_nonlinear_Data;

              redeclare function extends f_nonlinear  "Note that this function always sees the complete mass fraction vector" 
              protected
                MassFraction[nX] Xfull = if size(X, 1) == nX then X else cat(1, X, {1 - sum(X)});
                Real[nX] Y(unit = "mol/mol") = massToMoleFractions(if size(X, 1) == nX then X else cat(1, X, {1 - sum(X)}), data.MM) "Molar fractions";
              algorithm
                y := s_TX(x, Xfull) - sum(Xfull[i] * Modelica.Constants.R / MMX[i] * (if Xfull[i] < Modelica.Constants.eps then Y[i] else Modelica.Math.log(Y[i] * p / reference_p)) for i in 1:nX);
              end f_nonlinear;

              redeclare function extends solve  end solve;
            end Internal;
          algorithm
            T := Internal.solve(s, 200, 6000, p, Xfull, data[1]);
          end T_psX;
        end MixtureGasNasa;

        package FluidData  "Critical data, dipole moments and related data" 
          extends Modelica.Icons.Package;
          constant Modelica.Media.Interfaces.Types.IdealGas.FluidConstants N2(chemicalFormula = "N2", iupacName = "unknown", structureFormula = "unknown", casRegistryNumber = "7727-37-9", meltingPoint = 63.15, normalBoilingPoint = 77.35, criticalTemperature = 126.20, criticalPressure = 33.98e5, criticalMolarVolume = 90.10e-6, acentricFactor = 0.037, dipoleMoment = 0.0, molarMass = SingleGasesData.N2.MM, hasDipoleMoment = true, hasIdealGasHeatCapacity = true, hasCriticalData = true, hasAcentricFactor = true);
          constant Modelica.Media.Interfaces.Types.IdealGas.FluidConstants O2(chemicalFormula = "O2", iupacName = "unknown", structureFormula = "unknown", casRegistryNumber = "7782-44-7", meltingPoint = 54.36, normalBoilingPoint = 90.17, criticalTemperature = 154.58, criticalPressure = 50.43e5, criticalMolarVolume = 73.37e-6, acentricFactor = 0.022, dipoleMoment = 0.0, molarMass = SingleGasesData.O2.MM, hasDipoleMoment = true, hasIdealGasHeatCapacity = true, hasCriticalData = true, hasAcentricFactor = true);
          constant Modelica.Media.Interfaces.Types.IdealGas.FluidConstants CO2(chemicalFormula = "CO2", iupacName = "unknown", structureFormula = "unknown", casRegistryNumber = "124-38-9", meltingPoint = 216.58, normalBoilingPoint = -1.0, criticalTemperature = 304.12, criticalPressure = 73.74e5, criticalMolarVolume = 94.07e-6, acentricFactor = 0.225, dipoleMoment = 0.0, molarMass = SingleGasesData.CO2.MM, hasDipoleMoment = true, hasIdealGasHeatCapacity = true, hasCriticalData = true, hasAcentricFactor = true);
        end FluidData;

        package SingleGasesData  "Ideal gas data based on the NASA Glenn coefficients" 
          extends Modelica.Icons.Package;
          constant IdealGases.Common.DataRecord Ar(name = "Ar", MM = 0.039948, Hf = 0, H0 = 155137.3785921698, Tlimit = 1000, alow = {0, 0, 2.5, 0, 0, 0, 0}, blow = {-745.375, 4.37967491}, ahigh = {20.10538475, -0.05992661069999999, 2.500069401, -3.99214116e-008, 1.20527214e-011, -1.819015576e-015, 1.078576636e-019}, bhigh = {-744.993961, 4.37918011}, R = 208.1323720837088);
          constant IdealGases.Common.DataRecord CH4(name = "CH4", MM = 0.01604246, Hf = -4650159.63885838, H0 = 624355.7409524474, Tlimit = 1000, alow = {-176685.0998, 2786.18102, -12.0257785, 0.0391761929, -3.61905443e-005, 2.026853043e-008, -4.976705489999999e-012}, blow = {-23313.1436, 89.0432275}, ahigh = {3730042.76, -13835.01485, 20.49107091, -0.001961974759, 4.72731304e-007, -3.72881469e-011, 1.623737207e-015}, bhigh = {75320.6691, -121.9124889}, R = 518.2791167938085);
          constant IdealGases.Common.DataRecord CH3OH(name = "CH3OH", MM = 0.03204186, Hf = -6271171.523750494, H0 = 356885.5553329301, Tlimit = 1000, alow = {-241664.2886, 4032.14719, -20.46415436, 0.0690369807, -7.59893269e-005, 4.59820836e-008, -1.158706744e-011}, blow = {-44332.61169999999, 140.014219}, ahigh = {3411570.76, -13455.00201, 22.61407623, -0.002141029179, 3.73005054e-007, -3.49884639e-011, 1.366073444e-015}, bhigh = {56360.8156, -127.7814279}, R = 259.4878075117987);
          constant IdealGases.Common.DataRecord CO(name = "CO", MM = 0.0280101, Hf = -3946262.098314536, H0 = 309570.6191695138, Tlimit = 1000, alow = {14890.45326, -292.2285939, 5.72452717, -0.008176235030000001, 1.456903469e-005, -1.087746302e-008, 3.027941827e-012}, blow = {-13031.31878, -7.85924135}, ahigh = {461919.725, -1944.704863, 5.91671418, -0.0005664282830000001, 1.39881454e-007, -1.787680361e-011, 9.62093557e-016}, bhigh = {-2466.261084, -13.87413108}, R = 296.8383547363272);
          constant IdealGases.Common.DataRecord CO2(name = "CO2", MM = 0.0440095, Hf = -8941478.544405185, H0 = 212805.6215135368, Tlimit = 1000, alow = {49436.5054, -626.411601, 5.30172524, 0.002503813816, -2.127308728e-007, -7.68998878e-010, 2.849677801e-013}, blow = {-45281.9846, -7.04827944}, ahigh = {117696.2419, -1788.791477, 8.29152319, -9.22315678e-005, 4.86367688e-009, -1.891053312e-012, 6.330036589999999e-016}, bhigh = {-39083.5059, -26.52669281}, R = 188.9244822140674);
          constant IdealGases.Common.DataRecord C2H2_vinylidene(name = "C2H2_vinylidene", MM = 0.02603728, Hf = 15930556.80163212, H0 = 417638.4015534649, Tlimit = 1000, alow = {-14660.42239, 278.9475593, 1.276229776, 0.01395015463, -1.475702649e-005, 9.476298110000001e-009, -2.567602217e-012}, blow = {47361.1018, 16.58225704}, ahigh = {1940838.725, -6892.718150000001, 13.39582494, -0.0009368968669999999, 1.470804368e-007, -1.220040365e-011, 4.12239166e-016}, bhigh = {91071.1293, -63.3750293}, R = 319.3295152181795);
          constant IdealGases.Common.DataRecord C2H4(name = "C2H4", MM = 0.02805316, Hf = 1871446.924339362, H0 = 374955.5843263291, Tlimit = 1000, alow = {-116360.5836, 2554.85151, -16.09746428, 0.0662577932, -7.885081859999999e-005, 5.12522482e-008, -1.370340031e-011}, blow = {-6176.19107, 109.3338343}, ahigh = {3408763.67, -13748.47903, 23.65898074, -0.002423804419, 4.43139566e-007, -4.35268339e-011, 1.775410633e-015}, bhigh = {88204.2938, -137.1278108}, R = 296.3827247982046);
          constant IdealGases.Common.DataRecord C2H6(name = "C2H6", MM = 0.03006904, Hf = -2788633.890539904, H0 = 395476.3437741943, Tlimit = 1000, alow = {-186204.4161, 3406.19186, -19.51705092, 0.0756583559, -8.20417322e-005, 5.0611358e-008, -1.319281992e-011}, blow = {-27029.3289, 129.8140496}, ahigh = {5025782.13, -20330.22397, 33.2255293, -0.00383670341, 7.23840586e-007, -7.3191825e-011, 3.065468699e-015}, bhigh = {111596.395, -203.9410584}, R = 276.5127187299628);
          constant IdealGases.Common.DataRecord C2H5OH(name = "C2H5OH", MM = 0.04606844, Hf = -5100020.751733725, H0 = 315659.1801241805, Tlimit = 1000, alow = {-234279.1392, 4479.18055, -27.44817302, 0.1088679162, -0.0001305309334, 8.437346399999999e-008, -2.234559017e-011}, blow = {-50222.29, 176.4829211}, ahigh = {4694817.65, -19297.98213, 34.4758404, -0.00323616598, 5.78494772e-007, -5.56460027e-011, 2.2262264e-015}, bhigh = {86016.22709999999, -203.4801732}, R = 180.4808671619877);
          constant IdealGases.Common.DataRecord C3H6_propylene(name = "C3H6_propylene", MM = 0.04207974, Hf = 475288.1077687267, H0 = 322020.9535515191, Tlimit = 1000, alow = {-191246.2174, 3542.07424, -21.14878626, 0.0890148479, -0.0001001429154, 6.267959389999999e-008, -1.637870781e-011}, blow = {-15299.61824, 140.7641382}, ahigh = {5017620.34, -20860.84035, 36.4415634, -0.00388119117, 7.27867719e-007, -7.321204500000001e-011, 3.052176369e-015}, bhigh = {126124.5355, -219.5715757}, R = 197.588483198803);
          constant IdealGases.Common.DataRecord C3H8(name = "C3H8", MM = 0.04409562, Hf = -2373931.923397381, H0 = 334301.1845620949, Tlimit = 1000, alow = {-243314.4337, 4656.27081, -29.39466091, 0.1188952745, -0.0001376308269, 8.814823909999999e-008, -2.342987994e-011}, blow = {-35403.3527, 184.1749277}, ahigh = {6420731.680000001, -26597.91134, 45.3435684, -0.00502066392, 9.471216939999999e-007, -9.57540523e-011, 4.00967288e-015}, bhigh = {145558.2459, -281.8374734}, R = 188.5555073270316);
          constant IdealGases.Common.DataRecord C4H8_1_butene(name = "C4H8_1_butene", MM = 0.05610631999999999, Hf = -9624.584182316718, H0 = 305134.9651875226, Tlimit = 1000, alow = {-272149.2014, 5100.079250000001, -31.8378625, 0.1317754442, -0.0001527359339, 9.714761109999999e-008, -2.56020447e-011}, blow = {-25230.96386, 200.6932108}, ahigh = {6257948.609999999, -26603.76305, 47.6492005, -0.00438326711, 7.12883844e-007, -5.991020839999999e-011, 2.051753504e-015}, bhigh = {156925.2657, -291.3869761}, R = 148.1913623991023);
          constant IdealGases.Common.DataRecord C4H10_n_butane(name = "C4H10_n_butane", MM = 0.0581222, Hf = -2164233.28779709, H0 = 330832.0228759407, Tlimit = 1000, alow = {-317587.254, 6176.331819999999, -38.9156212, 0.1584654284, -0.0001860050159, 1.199676349e-007, -3.20167055e-011}, blow = {-45403.63390000001, 237.9488665}, ahigh = {7682322.45, -32560.5151, 57.3673275, -0.00619791681, 1.180186048e-006, -1.221893698e-010, 5.250635250000001e-015}, bhigh = {177452.656, -358.791876}, R = 143.0515706563069);
          constant IdealGases.Common.DataRecord C5H10_1_pentene(name = "C5H10_1_pentene", MM = 0.07013290000000001, Hf = -303423.9279995551, H0 = 309127.3852927798, Tlimit = 1000, alow = {-534054.813, 9298.917380000001, -56.6779245, 0.2123100266, -0.000257129829, 1.666834304e-007, -4.43408047e-011}, blow = {-47906.8218, 339.60364}, ahigh = {3744014.97, -21044.85321, 47.3612699, -0.00042442012, -3.89897505e-008, 1.367074243e-011, -9.31319423e-016}, bhigh = {115409.1373, -278.6177449000001}, R = 118.5530899192818);
          constant IdealGases.Common.DataRecord C5H12_n_pentane(name = "C5H12_n_pentane", MM = 0.07214878, Hf = -2034130.029641527, H0 = 335196.2430965569, Tlimit = 1000, alow = {-276889.4625, 5834.28347, -36.1754148, 0.1533339707, -0.0001528395882, 8.191092e-008, -1.792327902e-011}, blow = {-46653.7525, 226.5544053}, ahigh = {-2530779.286, -8972.59326, 45.3622326, -0.002626989916, 3.135136419e-006, -5.31872894e-010, 2.886896868e-014}, bhigh = {14846.16529, -251.6550384}, R = 115.2406457877736);
          constant IdealGases.Common.DataRecord C6H6(name = "C6H6", MM = 0.07811184, Hf = 1061042.730525872, H0 = 181735.4577743912, Tlimit = 1000, alow = {-167734.0902, 4404.50004, -37.1737791, 0.1640509559, -0.0002020812374, 1.307915264e-007, -3.4442841e-011}, blow = {-10354.55401, 216.9853345}, ahigh = {4538575.72, -22605.02547, 46.940073, -0.004206676830000001, 7.90799433e-007, -7.9683021e-011, 3.32821208e-015}, bhigh = {139146.4686, -286.8751333}, R = 106.4431717393932);
          constant IdealGases.Common.DataRecord C6H12_1_hexene(name = "C6H12_1_hexene", MM = 0.08415948000000001, Hf = -498458.4030224521, H0 = 311788.9986962847, Tlimit = 1000, alow = {-666883.165, 11768.64939, -72.70998330000001, 0.2709398396, -0.00033332464, 2.182347097e-007, -5.85946882e-011}, blow = {-62157.8054, 428.682564}, ahigh = {733290.696, -14488.48641, 46.7121549, 0.00317297847, -5.24264652e-007, 4.28035582e-011, -1.472353254e-015}, bhigh = {66977.4041, -262.3643854}, R = 98.79424159940152);
          constant IdealGases.Common.DataRecord C6H14_n_hexane(name = "C6H14_n_hexane", MM = 0.08617535999999999, Hf = -1936980.593988816, H0 = 333065.0431863586, Tlimit = 1000, alow = {-581592.67, 10790.97724, -66.3394703, 0.2523715155, -0.0002904344705, 1.802201514e-007, -4.617223680000001e-011}, blow = {-72715.4457, 393.828354}, ahigh = {-3106625.684, -7346.087920000001, 46.94131760000001, 0.001693963977, 2.068996667e-006, -4.21214168e-010, 2.452345845e-014}, bhigh = {523.750312, -254.9967718}, R = 96.48317105956971);
          constant IdealGases.Common.DataRecord C7H14_1_heptene(name = "C7H14_1_heptene", MM = 0.09818605999999999, Hf = -639194.6066478277, H0 = 313588.3036756949, Tlimit = 1000, alow = {-744940.284, 13321.79893, -82.81694379999999, 0.3108065994, -0.000378677992, 2.446841042e-007, -6.488763869999999e-011}, blow = {-72178.8501, 485.667149}, ahigh = {-1927608.174, -9125.024420000002, 47.4817797, 0.00606766053, -8.684859080000001e-007, 5.81399526e-011, -1.473979569e-015}, bhigh = {26009.14656, -256.2880707}, R = 84.68077851377274);
          constant IdealGases.Common.DataRecord C7H16_n_heptane(name = "C7H16_n_heptane", MM = 0.10020194, Hf = -1874015.612871368, H0 = 331540.487140269, Tlimit = 1000, alow = {-612743.289, 11840.85437, -74.87188599999999, 0.2918466052, -0.000341679549, 2.159285269e-007, -5.65585273e-011}, blow = {-80134.0894, 440.721332}, ahigh = {9135632.469999999, -39233.1969, 78.8978085, -0.00465425193, 2.071774142e-006, -3.4425393e-010, 1.976834775e-014}, bhigh = {205070.8295, -485.110402}, R = 82.97715593131233);
          constant IdealGases.Common.DataRecord C8H10_ethylbenz(name = "C8H10_ethylbenz", MM = 0.106165, Hf = 281825.4603682946, H0 = 209862.0072528611, Tlimit = 1000, alow = {-469494, 9307.16836, -65.2176947, 0.2612080237, -0.000318175348, 2.051355473e-007, -5.40181735e-011}, blow = {-40738.7021, 378.090436}, ahigh = {5551564.100000001, -28313.80598, 60.6124072, 0.001042112857, -1.327426719e-006, 2.166031743e-010, -1.142545514e-014}, bhigh = {164224.1062, -369.176982}, R = 78.31650732350586);
          constant IdealGases.Common.DataRecord C8H18_n_octane(name = "C8H18_n_octane", MM = 0.11422852, Hf = -1827477.060895125, H0 = 330740.51909278, Tlimit = 1000, alow = {-698664.715, 13385.01096, -84.1516592, 0.327193666, -0.000377720959, 2.339836988e-007, -6.01089265e-011}, blow = {-90262.2325, 493.922214}, ahigh = {6365406.949999999, -31053.64657, 69.6916234, 0.01048059637, -4.12962195e-006, 5.543226319999999e-010, -2.651436499e-014}, bhigh = {150096.8785, -416.989565}, R = 72.78805678301707);
          constant IdealGases.Common.DataRecord CL2(name = "CL2", MM = 0.07090600000000001, Hf = 0, H0 = 129482.8364313316, Tlimit = 1000, alow = {34628.1517, -554.7126520000001, 6.20758937, -0.002989632078, 3.17302729e-006, -1.793629562e-009, 4.260043590000001e-013}, blow = {1534.069331, -9.438331107}, ahigh = {6092569.42, -19496.27662, 28.54535795, -0.01449968764, 4.46389077e-006, -6.35852586e-010, 3.32736029e-014}, bhigh = {121211.7724, -169.0778824}, R = 117.2604857134798);
          constant IdealGases.Common.DataRecord F2(name = "F2", MM = 0.0379968064, Hf = 0, H0 = 232259.1511269747, Tlimit = 1000, alow = {10181.76308, 22.74241183, 1.97135304, 0.008151604010000001, -1.14896009e-005, 7.95865253e-009, -2.167079526e-012}, blow = {-958.6943, 11.30600296}, ahigh = {-2941167.79, 9456.5977, -7.73861615, 0.00764471299, -2.241007605e-006, 2.915845236e-010, -1.425033974e-014}, bhigh = {-60710.0561, 84.23835080000001}, R = 218.8202848542556);
          constant IdealGases.Common.DataRecord H2(name = "H2", MM = 0.00201588, Hf = 0, H0 = 4200697.462150524, Tlimit = 1000, alow = {40783.2321, -800.918604, 8.21470201, -0.01269714457, 1.753605076e-005, -1.20286027e-008, 3.36809349e-012}, blow = {2682.484665, -30.43788844}, ahigh = {560812.801, -837.150474, 2.975364532, 0.001252249124, -3.74071619e-007, 5.936625200000001e-011, -3.6069941e-015}, bhigh = {5339.82441, -2.202774769}, R = 4124.487568704486);
          constant IdealGases.Common.DataRecord H2O(name = "H2O", MM = 0.01801528, Hf = -13423382.81725291, H0 = 549760.6476280135, Tlimit = 1000, alow = {-39479.6083, 575.573102, 0.931782653, 0.00722271286, -7.34255737e-006, 4.95504349e-009, -1.336933246e-012}, blow = {-33039.7431, 17.24205775}, ahigh = {1034972.096, -2412.698562, 4.64611078, 0.002291998307, -6.836830479999999e-007, 9.426468930000001e-011, -4.82238053e-015}, bhigh = {-13842.86509, -7.97814851}, R = 461.5233290850878);
          constant IdealGases.Common.DataRecord He(name = "He", MM = 0.004002602, Hf = 0, H0 = 1548349.798456104, Tlimit = 1000, alow = {0, 0, 2.5, 0, 0, 0, 0}, blow = {-745.375, 0.9287239740000001}, ahigh = {0, 0, 2.5, 0, 0, 0, 0}, bhigh = {-745.375, 0.9287239740000001}, R = 2077.26673798694);
          constant IdealGases.Common.DataRecord NH3(name = "NH3", MM = 0.01703052, Hf = -2697510.117130892, H0 = 589713.1150428759, Tlimit = 1000, alow = {-76812.26149999999, 1270.951578, -3.89322913, 0.02145988418, -2.183766703e-005, 1.317385706e-008, -3.33232206e-012}, blow = {-12648.86413, 43.66014588}, ahigh = {2452389.535, -8040.89424, 12.71346201, -0.000398018658, 3.55250275e-008, 2.53092357e-012, -3.32270053e-016}, bhigh = {43861.91959999999, -64.62330602}, R = 488.2101075011215);
          constant IdealGases.Common.DataRecord NO(name = "NO", MM = 0.0300061, Hf = 3041758.509103149, H0 = 305908.1320131574, Tlimit = 1000, alow = {-11439.16503, 153.6467592, 3.43146873, -0.002668592368, 8.48139912e-006, -7.685111050000001e-009, 2.386797655e-012}, blow = {9098.214410000001, 6.72872549}, ahigh = {223901.8716, -1289.651623, 5.43393603, -0.00036560349, 9.880966450000001e-008, -1.416076856e-011, 9.380184619999999e-016}, bhigh = {17503.17656, -8.50166909}, R = 277.0927244793559);
          constant IdealGases.Common.DataRecord NO2(name = "NO2", MM = 0.0460055, Hf = 743237.6346306421, H0 = 221890.3174620426, Tlimit = 1000, alow = {-56420.3878, 963.308572, -2.434510974, 0.01927760886, -1.874559328e-005, 9.145497730000001e-009, -1.777647635e-012}, blow = {-1547.925037, 40.6785121}, ahigh = {721300.157, -3832.6152, 11.13963285, -0.002238062246, 6.54772343e-007, -7.6113359e-011, 3.32836105e-015}, bhigh = {25024.97403, -43.0513004}, R = 180.7277825477389);
          constant IdealGases.Common.DataRecord N2(name = "N2", MM = 0.0280134, Hf = 0, H0 = 309498.4543111511, Tlimit = 1000, alow = {22103.71497, -381.846182, 6.08273836, -0.00853091441, 1.384646189e-005, -9.62579362e-009, 2.519705809e-012}, blow = {710.846086, -10.76003744}, ahigh = {587712.406, -2239.249073, 6.06694922, -0.00061396855, 1.491806679e-007, -1.923105485e-011, 1.061954386e-015}, bhigh = {12832.10415, -15.86640027}, R = 296.8033869505308);
          constant IdealGases.Common.DataRecord N2O(name = "N2O", MM = 0.0440128, Hf = 1854006.107314236, H0 = 217685.1961247637, Tlimit = 1000, alow = {42882.2597, -644.011844, 6.03435143, 0.0002265394436, 3.47278285e-006, -3.62774864e-009, 1.137969552e-012}, blow = {11794.05506, -10.0312857}, ahigh = {343844.804, -2404.557558, 9.125636220000001, -0.000540166793, 1.315124031e-007, -1.4142151e-011, 6.38106687e-016}, bhigh = {21986.32638, -31.47805016}, R = 188.9103169986913);
          constant IdealGases.Common.DataRecord Ne(name = "Ne", MM = 0.0201797, Hf = 0, H0 = 307111.9986917546, Tlimit = 1000, alow = {0, 0, 2.5, 0, 0, 0, 0}, blow = {-745.375, 3.35532272}, ahigh = {0, 0, 2.5, 0, 0, 0, 0}, bhigh = {-745.375, 3.35532272}, R = 412.0215860493466);
          constant IdealGases.Common.DataRecord O2(name = "O2", MM = 0.0319988, Hf = 0, H0 = 271263.4223783392, Tlimit = 1000, alow = {-34255.6342, 484.700097, 1.119010961, 0.00429388924, -6.83630052e-007, -2.0233727e-009, 1.039040018e-012}, blow = {-3391.45487, 18.4969947}, ahigh = {-1037939.022, 2344.830282, 1.819732036, 0.001267847582, -2.188067988e-007, 2.053719572e-011, -8.193467050000001e-016}, bhigh = {-16890.10929, 17.38716506}, R = 259.8369938872708);
          constant IdealGases.Common.DataRecord SO2(name = "SO2", MM = 0.0640638, Hf = -4633037.690552231, H0 = 164650.3485587805, Tlimit = 1000, alow = {-53108.4214, 909.031167, -2.356891244, 0.02204449885, -2.510781471e-005, 1.446300484e-008, -3.36907094e-012}, blow = {-41137.52080000001, 40.45512519}, ahigh = {-112764.0116, -825.226138, 7.61617863, -0.000199932761, 5.65563143e-008, -5.45431661e-012, 2.918294102e-016}, bhigh = {-33513.0869, -16.55776085}, R = 129.7842463294403);
          constant IdealGases.Common.DataRecord SO3(name = "SO3", MM = 0.0800632, Hf = -4944843.573576874, H0 = 145990.9046852986, Tlimit = 1000, alow = {-39528.5529, 620.857257, -1.437731716, 0.02764126467, -3.144958662e-005, 1.792798e-008, -4.12638666e-012}, blow = {-51841.0617, 33.91331216}, ahigh = {-216692.3781, -1301.022399, 10.96287985, -0.000383710002, 8.466889039999999e-008, -9.70539929e-012, 4.49839754e-016}, bhigh = {-43982.83990000001, -36.55217314}, R = 103.8488594010732);
        end SingleGasesData;
      end Common;
    end IdealGases;
  end Media;

  package Math  "Library of mathematical functions (e.g., sin, cos) and of functions operating on vectors and matrices" 
    extends Modelica.Icons.Package;

    package Icons  "Icons for Math" 
      extends Modelica.Icons.IconsPackage;

      partial function AxisLeft  "Basic icon for mathematical function with y-axis on left side" end AxisLeft;

      partial function AxisCenter  "Basic icon for mathematical function with y-axis in the center" end AxisCenter;
    end Icons;

    function asin  "Inverse sine (-1 <= u <= 1)" 
      extends Modelica.Math.Icons.AxisCenter;
      input Real u;
      output .Modelica.SIunits.Angle y;
      external "builtin" y = asin(u);
    end asin;

    function exp  "Exponential, base e" 
      extends Modelica.Math.Icons.AxisCenter;
      input Real u;
      output Real y;
      external "builtin" y = exp(u);
    end exp;

    function log  "Natural (base e) logarithm (u shall be > 0)" 
      extends Modelica.Math.Icons.AxisLeft;
      input Real u;
      output Real y;
      external "builtin" y = log(u);
    end log;
  end Math;

  package Utilities  "Library of utility functions dedicated to scripting (operating on files, streams, strings, system)" 
    extends Modelica.Icons.Package;

    package Streams  "Read from files and write to files" 
      extends Modelica.Icons.Package;

      function error  "Print error message and cancel all actions" 
        extends Modelica.Icons.Function;
        input String string "String to be printed to error message window";
        external "C" ModelicaError(string) annotation(Library = "ModelicaExternalC");
      end error;
    end Streams;
  end Utilities;

  package Constants  "Library of mathematical constants and constants of nature (e.g., pi, eps, R, sigma)" 
    extends Modelica.Icons.Package;
    final constant Real pi = 2 * Math.asin(1.0);
    final constant Real eps = ModelicaServices.Machine.eps "Biggest number such that 1.0 + eps = 1.0";
    final constant .Modelica.SIunits.Velocity c = 299792458 "Speed of light in vacuum";
    final constant Real R(final unit = "J/(mol.K)") = 8.3144598 "Molar gas constant (previous value: 8.314472)";
    final constant Real mue_0(final unit = "N/A2") = 4 * pi * 1.e-7 "Magnetic constant";
    final constant .Modelica.SIunits.Conversions.NonSIunits.Temperature_degC T_zero = -273.15 "Absolute zero temperature";
  end Constants;

  package Icons  "Library of icons" 
    extends Icons.Package;

    partial package Package  "Icon for standard packages" end Package;

    partial package VariantsPackage  "Icon for package containing variants" 
      extends Modelica.Icons.Package;
    end VariantsPackage;

    partial package InterfacesPackage  "Icon for packages containing interfaces" 
      extends Modelica.Icons.Package;
    end InterfacesPackage;

    partial package UtilitiesPackage  "Icon for utility packages" 
      extends Modelica.Icons.Package;
    end UtilitiesPackage;

    partial package IconsPackage  "Icon for packages containing icons" 
      extends Modelica.Icons.Package;
    end IconsPackage;

    partial package MaterialPropertiesPackage  "Icon for package containing property classes" 
      extends Modelica.Icons.Package;
    end MaterialPropertiesPackage;

    partial function Function  "Icon for functions" end Function;

    partial record Record  "Icon for records" end Record;
  end Icons;

  package SIunits  "Library of type and unit definitions based on SI units according to ISO 31-1992" 
    extends Modelica.Icons.Package;

    package Icons  "Icons for SIunits" 
      extends Modelica.Icons.IconsPackage;

      partial function Conversion  "Base icon for conversion functions" end Conversion;
    end Icons;

    package Conversions  "Conversion functions to/from non SI units and type definitions of non SI units" 
      extends Modelica.Icons.Package;

      package NonSIunits  "Type definitions of non SI units" 
        extends Modelica.Icons.Package;
        type Temperature_degC = Real(final quantity = "ThermodynamicTemperature", final unit = "degC") "Absolute temperature in degree Celsius (for relative temperature use SIunits.TemperatureDifference)" annotation(absoluteValue = true);
        type Pressure_bar = Real(final quantity = "Pressure", final unit = "bar") "Absolute pressure in bar";
      end NonSIunits;

      function to_degC  "Convert from Kelvin to degCelsius" 
        extends Modelica.SIunits.Icons.Conversion;
        input Temperature Kelvin "Kelvin value";
        output NonSIunits.Temperature_degC Celsius "Celsius value";
      algorithm
        Celsius := Kelvin + Modelica.Constants.T_zero;
        annotation(Inline = true); 
      end to_degC;

      function from_degC  "Convert from degCelsius to Kelvin" 
        extends Modelica.SIunits.Icons.Conversion;
        input NonSIunits.Temperature_degC Celsius "Celsius value";
        output Temperature Kelvin "Kelvin value";
      algorithm
        Kelvin := Celsius - Modelica.Constants.T_zero;
        annotation(Inline = true); 
      end from_degC;

      function to_bar  "Convert from Pascal to bar" 
        extends Modelica.SIunits.Icons.Conversion;
        input Pressure Pa "Pascal value";
        output NonSIunits.Pressure_bar bar "bar value";
      algorithm
        bar := Pa / 1e5;
        annotation(Inline = true); 
      end to_bar;
    end Conversions;

    type Angle = Real(final quantity = "Angle", final unit = "rad", displayUnit = "deg");
    type Area = Real(final quantity = "Area", final unit = "m2");
    type Volume = Real(final quantity = "Volume", final unit = "m3");
    type Velocity = Real(final quantity = "Velocity", final unit = "m/s");
    type Acceleration = Real(final quantity = "Acceleration", final unit = "m/s2");
    type Mass = Real(quantity = "Mass", final unit = "kg", min = 0);
    type Density = Real(final quantity = "Density", final unit = "kg/m3", displayUnit = "g/cm3", min = 0.0);
    type Pressure = Real(final quantity = "Pressure", final unit = "Pa", displayUnit = "bar");
    type AbsolutePressure = Pressure(min = 0.0, nominal = 1e5);
    type DynamicViscosity = Real(final quantity = "DynamicViscosity", final unit = "Pa.s", min = 0);
    type Energy = Real(final quantity = "Energy", final unit = "J");
    type Power = Real(final quantity = "Power", final unit = "W");
    type MassFlowRate = Real(quantity = "MassFlowRate", final unit = "kg/s");
    type MomentumFlux = Real(final quantity = "MomentumFlux", final unit = "N");
    type ThermodynamicTemperature = Real(final quantity = "ThermodynamicTemperature", final unit = "K", min = 0.0, start = 288.15, nominal = 300, displayUnit = "degC") "Absolute temperature (use type TemperatureDifference for relative temperatures)" annotation(absoluteValue = true);
    type Temp_K = ThermodynamicTemperature;
    type Temperature = ThermodynamicTemperature;
    type Compressibility = Real(final quantity = "Compressibility", final unit = "1/Pa");
    type IsothermalCompressibility = Compressibility;
    type ThermalConductivity = Real(final quantity = "ThermalConductivity", final unit = "W/(m.K)");
    type SpecificHeatCapacity = Real(final quantity = "SpecificHeatCapacity", final unit = "J/(kg.K)");
    type RatioOfSpecificHeatCapacities = Real(final quantity = "RatioOfSpecificHeatCapacities", final unit = "1");
    type Entropy = Real(final quantity = "Entropy", final unit = "J/K");
    type SpecificEntropy = Real(final quantity = "SpecificEntropy", final unit = "J/(kg.K)");
    type SpecificEnergy = Real(final quantity = "SpecificEnergy", final unit = "J/kg");
    type SpecificEnthalpy = SpecificEnergy;
    type DerDensityByPressure = Real(final unit = "s2/m2");
    type DerDensityByTemperature = Real(final unit = "kg/(m3.K)");
    type AmountOfSubstance = Real(final quantity = "AmountOfSubstance", final unit = "mol", min = 0);
    type MolarMass = Real(final quantity = "MolarMass", final unit = "kg/mol", min = 0);
    type MolarVolume = Real(final quantity = "MolarVolume", final unit = "m3/mol", min = 0);
    type MassFraction = Real(final quantity = "MassFraction", final unit = "1", min = 0, max = 1);
    type MoleFraction = Real(final quantity = "MoleFraction", final unit = "1", min = 0, max = 1);
    type FaradayConstant = Real(final quantity = "FaradayConstant", final unit = "C/mol");
  end SIunits;
  annotation(version = "3.2.2", versionBuild = 3, versionDate = "2016-04-03", dateModified = "2016-04-03 08:44:41Z"); 
end Modelica;

model plant_total
  extends plant;
end plant_total;
