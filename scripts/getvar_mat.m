function var = getvar_mat(varname)

load work/plant_res.mat ;
names = name';
dat   = data_2' ;

indx = -1 ;
for i = 1 : size(names,1)
    if (names(i, 1:length(varname)) == varname)
        indx = i;
    end

if (indx != -1) 
    var = dat(:, indx);
else
    var = [] ;
end

end 
