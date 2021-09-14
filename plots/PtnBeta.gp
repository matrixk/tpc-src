dfname="PtnBeta67.dat"
snoise=30.0

set output "1.eps"
set terminal postscript eps enhanced color size 6.25,2.3 "Helvetica" 12

#set cbrange [-3.0:10.0]

## CERN ROOT palette style 1
#set palette model RGB defined (0.0 1.0 1.0 1.0, 1/256. 1.0 1.0 1.0, 1/255. 0.0 0.0 0.51, 0.34 0.0 0.81 1.0, 0.61 0.87 1.0 0.12, 0.84 1.0 0.2 0.0, 1.0 0.51 0.0 0.0) positive

set size ratio -1
set view equal xy
set xrange [-70:70]
set yrange [-70:70]
unset zeroaxis
unset xlabel
unset ylabel
unset colorbox
#set format x ""
set cbrange [-300:3000]
set cblabel 'Number of charges per pixel'
set border lw 1.0

set multiplot layout 1,3 title "" \
    margins screen 0.06, 0.92, 0.15, 0.9 spacing screen 0.0,0.0

set xlabel 'x [mm]'
set ylabel 'y [mm]'

set title 'Pixel pitch = 1 [mm]'

plot dfname index 0:9918 u 6:7:($8+snoise*invnorm(rand(0))) w filledcurves fillcolor palette z t ''

unset xlabel
unset ylabel
set format y ""

set title '4 [mm]'

plot dfname index 10000:10720 u 6:7:($8+snoise*invnorm(rand(0))) w filledcurves fillcolor palette z t ''

set title '8 [mm]'
set colorbox
plot dfname index 20000:20216 u 6:7:($8+snoise*invnorm(rand(0))) w filledcurves fillcolor palette z t ''

# unset title
# unset colorbox
# snoise=0.0
# set xlabel 'x [mm]'
# set format x "%g"

# plot dfname index 0:9918 u 6:7:($8+snoise*invnorm(rand(0))) w filledcurves fillcolor palette z t ''

# unset xlabel
# unset ylabel
# set format y ""

# plot dfname index 10000:10720 u 6:7:($8+snoise*invnorm(rand(0))) w filledcurves fillcolor palette z t ''

# plot dfname index 20000:20216 u 6:7:($8+snoise*invnorm(rand(0))) w filledcurves fillcolor palette z t ''

unset multiplot

unset output
