normal(x, mu, sd) = (1/(sd*sqrt(2*pi)))*exp(-(x-mu)**2/(2*sd**2))
stats "data/test.dat"  u 1 every ::10000
bin_width =  0.2
set boxwidth bin_width
bin_number(x) = floor(x/bin_width)
rounded(x) = bin_width * ( bin_number(x) + 0.5 )
set terminal postscript color
set output "plots/test.ps"

plot 	"data/test.dat" every ::10000 using (rounded($1)):(1/(bin_width*STATS_records)) smooth frequency with boxes lc rgb "red" notitle,\
	(normal(x,0,1)+normal(x,10,2))/2
