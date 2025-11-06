#Some functions
#/degC(Rankine)=(Rankine - 491.67)/1.8 
#/bar(psi) = psi / 14.7

dataFile = "./work/plant_res.csv"
set datafile separator ","
set key autotitle columnhead

set xlabel "time (s)"
set multiplot layout 2,1

plot \
 dataFile   using (column("time")):(column("v1.inlet.m_flow")*60e3) \
             smooth cspline lw 3, \
 dataFile   using (column("time")):(column("v2.inlet.m_flow")*60e3) \
             smooth cspline lw 3 
set ylabel "Flow Rate (lpm)"
#set ytics 0.4
set ylabel "Valve Openining %"
plot  \
 dataFile   using 1:(column("v1.po")) \
             smooth cspline  lw 3, \
 dataFile   using 1:(column("v2.po")) \
             smooth cspline  lw 3 
unset object 1
unset multiplot
pause -1 
#
