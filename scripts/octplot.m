setenv("GNUTERM","wxt");
set(0, "defaultlinelinewidth", 2);
clear;clc;close all;
h = figure(1) ;
set(h,'resizefcn','refresh');
FS = findall(h, '-property', 'FontSize');
set(FS,'FontSize',6);

% Import Simulation Data
omimport("work/plant");

subplot(2,2,1);
plot(time, cmp.pratio,";pratio;");
grid on ;
xlabel("time");

subplot(2,2,2);
plot(time, cmp.speed, ";speed;");
grid on ;
xlabel("time");

subplot(2,2,3);
plot( time, cmp.state_out.T, ";cmp_T_out;",
      time, heater.Tf, ";heater_Tf;",
      time, heater.Tw, ";heater_Tw;");
grid on ;
xlabel("time");

subplot(2,2,4);
plot(time, hdrop, time, vdrop)
legend("heater_dp", "valve_dp");
grid on ;
xlabel("time");

pause;
