dfname="dpn_n1000.dat"

set output "1.eps"
set terminal postscript eps enhanced color size 6.25,2.0 "Helvetica" 12

#set cbrange [-3.0:10.0]

## CERN ROOT palette style 1
#set palette model RGB defined (0.0 1.0 1.0 1.0, 1/256. 1.0 1.0 1.0, 1/255. 0.0 0.0 0.51, 0.34 0.0 0.81 1.0, 0.61 0.87 1.0 0.12, 0.84 1.0 0.2 0.0, 1.0 0.51 0.0 0.0) positive

#set size ratio -1
#set view equal xy
set xrange [0:11]
set yrange [0:11]
unset xlabel
unset ylabel
unset colorbox
set view equal xy
set border lw 1.0

set multiplot layout 1,3 title "" \
    margins screen 0.08, 0.90, 0.15, 0.9 spacing screen 0.01,0.02

set xlabel '{/Symbol s}_t [mm]'
set ylabel 'Pixel pitch [mm]'
set cblabel 'Energy resolution FWHM [%]'
set cbrange [0:3]

set title 'ENC = 10 e^-/pixel'

plot sprintf("<awk '{if($3==10.0)print}' %s", dfname) u 1:2:($9/$8*2.35*100) w p pt 5 ps 2 lc palette t ''

unset xlabel
unset ylabel
set format y ""

set title '20 e^-/pixel'

plot sprintf("<awk '{if($3==20.0)print}' %s", dfname) u 1:2:($9/$8*2.35*100) w p pt 5 ps 2 lc palette t ''

set title '30 e^-/pixel'
set colorbox
plot sprintf("<awk '{if($3==30.0)print}' %s", dfname) u 1:2:($9/$8*2.35*100) w p pt 5 ps 2 lc palette t ''

# unset title
# unset colorbox
# set cbrange [0:200]
# set cblabel "Number of pixels with a hit"
# set format x "%g"
# set xlabel '{/Symbol s}_t [mm]'

# set palette model RGB defined (0.0 1.0 1.0 1.0, 1/256. 1.0 1.0 1.0, 1/255. 0.0 0.0 0.51, 0.34 0.0 0.81 1.0, 0.61 0.87 1.0 0.12, 0.84 1.0 0.2 0.0, 1.0 0.51 0.0 0.0) positive

# set format y "%g"
# plot sprintf("<awk '{if($3==10.0)print}' %s", dfname) u 1:2:6 w p pt 5 ps 2 lc palette t ''

# unset xlabel
# set format y ""
# plot sprintf("<awk '{if($3==20.0)print}' %s", dfname) u 1:2:6 w p pt 5 ps 2 lc palette t ''

# set colorbox

# plot sprintf("<awk '{if($3==30.0)print}' %s", dfname) u 1:2:6 w p pt 5 ps 2 lc palette t ''

unset multiplot

unset output
