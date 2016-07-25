function foo
input Real vec[:] = {1, 2, 3} ;
output Integer nV ;
    
algorithm
        nV := size(vec,1) ; 
end foo;
