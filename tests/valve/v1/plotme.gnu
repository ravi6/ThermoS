#Some functions
#/degC(Rankine)=(Rankine - 491.67)/1.8 
#/bar(psi) = psi / 14.7

dataFile = "./work/plant_res.csv"
set datafile separator ","
set key autotitle columnhead

set xlabel "time (s)"
set multiplot layout 2,1

plot for [i = 1:3] \
 dataFile   using 1:(column(sprintf("valve[%d].inlet.m_flow",i))) \
             smooth cspline ls i  lw 3 t sprintf("valve %d", i) 
set ylabel "Flow Rate (kg/s)"
#set ytics 0.4
set ylabel "Valve Openining %"

plot for [i = 1:3] \
 dataFile   using 1:(column(sprintf("valve[%d].po",i))) \
             smooth cspline ls i lw 3 t ""
unset object 1
unset multiplot
pause -1 

#
