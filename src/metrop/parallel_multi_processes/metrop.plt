#===================================================
# metrop.plt
#
# @Author: Wenchong Chen
#
# It generates figures
# - X = <m^2>
# - autocorrelation function Rho(t)
# - integrated autocorrelation time Tau(t)
# for J/K(B)T = log(1 + sqrt(2)) / 2
#===================================================


# setting shared by all figures
set terminal png

#============= X = <m^2> =============#
reset

n=32      					# number of bins
max=1.	  					# max value on x-axis
min=0.	  					# min value on x-axis
width=(max-min)/n  	# width of a bin

# map values to the intervals
hist(x,width)=width*floor(x/width)+width/2.0

set xrange [min:max]
set yrange [0:1400]

set xtics min,(max-min)/5,max
set style fill solid 0.5

set xlabel "X = <m^2>"
set ylabel "Frequency"

#----- Susceptibility with J/K(B)T = log(1 + sqrt(2)) / 2 -----#
set title "Susceptibility with J/K(B)T = log(1 + sqrt(2)) / 2"

set output "xt.png"

plot "xt.dat" u (hist($1,width)):(1.0) smooth freq w boxes lc rgb"green" notitle

unset output


#======= Autocorrelation Function Rho(t) =======#
reset

set key left

set xlabel "Time"
set ylabel "Rho(t)"

#----- Rho(t) with J/K(B)T = log(1 + sqrt(2)) / 2 -----#
set title "Autocorrelation Function with J/K(B)T = log(1 + sqrt(2)) / 2"

set output "rho.png"

plot "rho.dat" with lines lc rgb "red"

unset output


#======= Integrated Autocorrelation Time Tau(t) =======#
reset

set key left

set xlabel "Time"
set ylabel "Tau(t)"

#----- Tau(t) with J/K(B)T = log(1 + sqrt(2)) / 2 -----#
set title "Integrated Autocorrelation Time with J/K(B)T = log(1 + sqrt(2)) / 2"

set output "tau.png"

plot "tau.dat" with lines lc rgb "red"

unset output

