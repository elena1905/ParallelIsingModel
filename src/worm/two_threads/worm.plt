#===================================================
# metrop.plt
#
# @Author: Wenchong Chen
#
# It generates figures
# - L^(7/4) / X
# - autocorrelation function Rho(t)
# - integrated autocorrelation time Tau(t)
# for J/K(B)T = log(1 + sqrt(2)) / 2
#===================================================


# setting shared by all figures
set terminal png

#============= L^(7/4) / X =============#
reset

n=35      					# number of bins
max=1.	  					# max value on x-axis
min=0.	  					# min value on x-axis
width=(max-min)/n  	# width of a bin

# map values to the intervals
hist(x,width)=width*floor(x/width)+width/2.0

set xrange [min:max]
set yrange [0:]

set xtics min,(max-min)/5,max
set style fill solid 0.5

set xlabel "L^(7/4) / X = (nKick/N) * L^(7/4)"
set ylabel "Frequency"

#----- L^(7/4) / X with J/K(B)T = log(1 + sqrt(2)) / 2 -----#
set title "L^(7/4) / X with J/K(B)T = log(1 + sqrt(2)) / 2"

set output "xt.png"

plot "xt.dat" u (hist($1,width)) smooth freq w boxes lc rgb"green" notitle

unset output


#======= Autocorrelation Function Rho(t) =======#
reset

set key left

set xlabel "Time"
set ylabel "Rho(t)"

#----- Rho(t) with J/K(B)T = log(1 + sqrt(2)) / 2 -----#
set title "Autocorrelation Function Rho(t) with J/K(B)T = log(1 + sqrt(2)) / 2"

set output "rho.png"

plot "rho.dat" with lines lc rgb "red"

unset output


#======= Integrated Autocorrelation Time Tau(t) =======#
reset

set key left

set xlabel "Time"
set ylabel "Tau(t)"

#----- Tau(t) with J/K(B)T = log(1 + sqrt(2)) / 2 -----#
set title "Integrated Autocorrelation Time Tau(t) with J/K(B)T = log(1 + sqrt(2)) / 2"

set output "tau.png"

plot "tau.dat" with lines lc rgb "red"

unset output

