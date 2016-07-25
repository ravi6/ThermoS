model plant
   Integer n ; 
   constant Integer subvec[4] = {2, 10, 20, 30};
   constant Integer svals[3,4] = { subvec  for k in 1:3 } ;
   Integer  z[3,4](each start=svals);
equation
    n = foo({1, 4});     
    assert(false,"Vector Size that is passed = " + String(n),AssertionLevel.warning);
end plant;
