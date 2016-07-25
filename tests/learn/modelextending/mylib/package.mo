within ;

package mylib

    model vehicle  "A vehicle"
       parameter Real speed = 10 ;
       parameter Real k = 2 ;
       Real x ;
       equation
       der(x) = speed ;
    end vehicle ;

    model bus  
       extends vehicle  ;
       Real y ;

       equation
        der(y) = speed * k ;
    end bus;

end mylib;
