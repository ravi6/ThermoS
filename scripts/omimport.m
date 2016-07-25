function omimport(modelname)  
%
% Read OpenModelica Result File into Workspace
%
%
% Feedback/problems: Christian Schaad, ingenieurbuero@christian-schaad.de
% Modified by R. Saripalli, 
%            to handle 
%	     Array variables of form X[i] with  (renamed to X_i)  
%            and to discard trailing non-ascii chars from variable names.
%  Tested on Octave  
%  4th Dec. 2013


load ([modelname,'_res.mat']);

%Sort out double times 
size(data_2)
size(dataInfo) 

times = data_2(1,:); 
%size(times)
deltat0=find(diff(times)==0) %Vector of indecides with duplicate time stamp
disp(['Duplicate values: ',num2str(length(deltat0)),'/',num2str(length(data_2(1,:)))])
%times = unique(times);

assignin('base','data_2',data_2);
assignin('base','dataInfo',dataInfo);
assignin('base','name',name);
assignin('base','deltat0',deltat0);

name=name';
size(name)
disp(['Number of Variables:', num2str(length(name))]);
for i=1:length(name)
        tmp =  strfind( name(i,:),  'der('  );
	if (isempty(tmp))    % when der is absent

                [S,iendchar]=(regexp(name(i,:),"^[-_.a-zA-Z0-9 [\\]]*"));
                namex = [name(i,:)] 

		  %dataInfo second row holds row index of data_2
                  %  which holds variable i data. 
		  %  its sign is assigned to data values. (why all this fuss??)
		  %  and this data is assigned to varVals.

                try
                     
	        	if dataInfo(2,i) <0  %dataInfo second row holds row index of data_2
	           	   assignin('base','varVals',-data_2(-dataInfo(2,i),:));
	                else
	                   assignin( 'base', 'varVals', data_2(dataInfo(2,i),:) );
                        endif  
                catch
                        {i, name(i,:), "Value assignment error ", dataInfo(2,i)}
                end_try_catch

	%	evalin( 'base', (['varVals(deltat0)=[];']) )  %We ignore the data if dup. time stamp

		varName = num2str(name(i,1:iendchar));  	% Grab the variable name
		varName = strrep(varName, "[","_"); 
                varName = strrep(varName, "]","");  %Handle variables of array type
		%varName = strrep(varName, "[","("); 
                %varName = strrep(varName, "]",",:)");  %Handle variables of array type

                try
		    evalin( 'base',[varName,'=varVals;'] );
                catch
                    ["problem in evalin >>",  varName]
             end_try_catch
	endif 
endfor  % end of loop over all variables

clear data_1 data_2 Aclass description modelname i dataInfo varVals deltat0 varName;
evalin('base',(['clear name data_2 dataInfo iendchar varVals']));

endfunction
