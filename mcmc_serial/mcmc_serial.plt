normal(x, mu, sd) = (1/(sd*sqrt(2*pi)))*exp(-(x-mu)**2/(2*sd**2))
set xrange [-5:20]
set title "Probability distribution P(x)"
set grid
set term postscript color #output postscript
set output "plots/mcmc_serial.ps"

n=100 #number of intervals
max=5. #max value
min=-5. #min value
width=(max-min)/n #interval width
#function used to map a value to the intervals
hist(x,width)=width*floor(x/width)+width/2.0
set boxwidth width*0.9
set style fill solid 0.5 #fillstyleset xlabel "x"
set ylabel "Frequency"
set xlabel "x"
plot "./data/test.dat" u (hist($1,width)):(1.0) smooth freq w boxes lc rgb"green" notitle
