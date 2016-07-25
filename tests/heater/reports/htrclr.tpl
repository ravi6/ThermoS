function show_htrclr
algorithm
OpenModelica.Scripting.print(
     "Summary:  uop "  + "\n" + 
     "Mc = " +        ThermoS.Util.decimalS(val(uop.Mc, tSample),3)                + "\t" + 
     "Tf = " +        ThermoS.Util.decimalS(val(uop.Tf, tSample),0)                + "\n" +
     "Q_ew = " +      ThermoS.Util.decimalS(val(uop.Q_ew, tSample),0)              + "\t" +
     "Q_wf = " +      ThermoS.Util.decimalS(val(uop.Q_wf, tSample),0)              + "\n" +
     "Pdrop = " +     ThermoS.Util.decimalS((val(uop.inlet.p, tSample)
                                -val(uop.outlet.p, tSample))*1e-5, 2)   + "\t" + 
     "Tw = " +        ThermoS.Util.decimalS(val(uop.Tw, tSample),0)                + "\t" +
     "Cp = " +        ThermoS.Util.decimalS(val(uop.Cp, tSample),0)                + "\n" 
     );
end show_htrclr;
