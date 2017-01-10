load "./work/plant_res.mat";
 m1=data_2(12,:);
 m2=data_2(13,:);
 m3=data_2(14,:);
 v=data_2(20,:);
 plot(v,m1)
 plot(v,[m1; m2; m3])
 xlabel("Valve Openning %")
 xlabel("% Open")
 ylabel("Mass Flow Rate (kg/s)")
 legend("Linear", "Equi-percent", "Quick-Opening")
pause;
