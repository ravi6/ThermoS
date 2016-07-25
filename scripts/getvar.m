function var = getvar(simData, varname)

# Get the index of the variable you are after
j = 0 ;
for i = 1:simData.nVars
   str = cstrcat("\"",varname,"\"");
   if strcmp(str, simData.names(1,i))
      j = i ; 
   end
end
 if j == 0 
     printf ("Variable: %s not found", varname);
     exit ;
 end
 var = simData.data(:,j) ;
end 
