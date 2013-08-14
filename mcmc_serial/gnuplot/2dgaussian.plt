normal(x,y,mux,sdx,muy,sdy) = (1/(sdx*sdy*sqrt(2*pi)))*exp(-(x-mux)**2/(2*sdx**2))*exp(-(y-muy)**2/(2*sdx**2))
set xrange [-30:30]
set yrange [-30:30]
set isosamples 200 
set title "Probabilty Distribution P(x,y)"
set view map
unset surface
#set cntrparam levels 0
set contour base
#set key outside
set pm3d
#set grid
set term postscript color
set output "gaussian2d.ps"
func(x,y) = normal(x,y,-3,6,5,6)+normal(x,y,-10,5,-10,5) 
splot func(x,y) title "" 

