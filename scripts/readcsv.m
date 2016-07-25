function simData =  readcsv(fname)

# Reads Modelica Simulation data in CSV format
#   returns a Data structure that is easy to process
#   in Octave or Matlab

# Sanitize Variable Names that are Two dimensional Vectors
#   The sed command we want to execute is
#     sed -r '1,1s/(\[[0-9]*),/\1_/g   fname  > /tmp/tmpfile
#     Well we need to add some more escapes b'cause we are passing through system

system("rm -f /tmp/tmpfile");
cmd = cstrcat("sed -r '1,1s/(\\[[0-9]*),/\\1_/g' " ,  fname ,  " > /tmp/tmpfile");

if  system(cmd) == 0   
    cmd="Success";
    simData.data = csvread(fname);
else
    printf("Error parsing %s \n", fname) ; exit ;
end

# Now parse the variable names

fid = fopen("/tmp/tmpfile");
stNames=fgetl(fid) ;  
fclose(fid);
simData.names=ostrsplit(stNames,',');
simData.nVars=size(simData.names,2);

end 
