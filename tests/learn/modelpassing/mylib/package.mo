within ;

package mylib

    model vehicle  "A vehicle"
       parameter Real speed = 10 ;
       Real x ;
       equation
       der(x) = speed ;
    end vehicle ;

    model bus  
       extends vehicle  ;
       Real y ;
       equation
        der(y) = speed ;
    end bus;

    package p  "some package"
        function doit
           output Real z ;
           algorithm
               z:=55;
        end doit;
    end p;

    

    model geep
   //     replaceable package v=p ;
       replaceable  model  v = vehicle(speed=4); 
       v myv;  // instantiate vehicle
       Real z ;
       equation
         der(z) = myv.x ;  // we should get quadratic output
                           //  myv.x is linear in t
    end geep;
end mylib;
