# This is a good way to make default line styles you want
# define line styles using explicit rgbcolor names
# (you can add point style in it as well)
# Here I chose filled circles as my markder type
# (I prefer colors to distinguish different data sets)

#Some functions
degC(Rankine)=(Rankine - 491.67)/1.8 
bar(psi) = psi / 14.7


#set output "my.ps"
set style line 1 lt 2 lc rgb "red" lw 1     ps 1 pt 6 
set style line 2 lt 2 lc rgb "orange" lw 1  ps 1 pt 6
set style line 3 lt 2 lc rgb "green" lw 1   ps 1 pt 6
set style line 5 lt 2 lc rgb "magenta" lw 1  ps 1 pt 6
set terminal wxt size 800,600 enhanced font 'Arial,8' persist
#set terminal postscript 

# Annotate and Beautify
# All plots against speed
#
# Use these when you need time formatted X scale
#set timefmt "%H:%MM:%SS"
#set xdata time
#set format x timefmt 
set xtics rotate

#set term wxt 1 
set grid
set key inside top left
set xlabel "speed (rad/s)"

set multiplot layout 2,2

set ylabel "Heat Input (kW)"
set ytics 40
plot "ssdata"   using 1:2  smooth cspline ls 1 t "", \
     "ssdata"   using 1:2  with p ls 1 t "" 

set ylabel "Pressure Ratio"
set ytics 0.4
plot "ssdata"   using 1:3  smooth cspline ls 1 t "", \
     "ssdata"   using 1:3  with p ls 1 t "" 

set ylabel "Mass Flow Rate (kg/s)"
set ytics 0.1
plot "ssdata"   using 1:4  smooth cspline ls 1 t "", \
     "ssdata"   using 1:4  with p ls 1 t ""  

set ylabel "Beta Values"
set ytics 0.2
plot "ssdata"   using 1:5  t "Comp" smooth cspline ls 1, \
     "ssdata"   using 1:5  with p ls 1 t "" , \
     "ssdata"   using 1:6  t "Turb" smooth cspline ls 2, \
     "ssdata"   using 1:6  with p ls 2 t ""
unset multiplot
pause -1 

#pause -1
#
