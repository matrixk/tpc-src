dfname="dpn_n1000.dat"

set output "1.eps"
set terminal postscript eps enhanced color size 6.0,2.0 "Helvetica" 12

#set cbrange [-3.0:10.0]

## CERN ROOT palette style 1
#set palette model RGB defined (0.0 1.0 1.0 1.0, 1/256. 1.0 1.0 1.0, 1/255. 0.0 0.0 0.51, 0.34 0.0 0.81 1.0, 0.61 0.87 1.0 0.12, 0.84 1.0 0.2 0.0, 1.0 0.51 0.0 0.0) positive

#set size ratio -1
#set view equal xy
#set view equal xy
set border lw 1.0

set multiplot layout 1,2 title ""
#    margins screen 0.08, 0.95, 0.15, 0.9 spacing screen 0.2,0.02

set origin 0,0
set size 0.55,1.0

set xrange [0:600]
set xlabel "Number of pixels with a hit"
set yrange [0:2.5]
set ytics nomirror
set ylabel "Energy resolution FWHM [%]"
set y2range [500:11000]
set y2tics
set my2tics 10
set format y2 "%g"
set y2label "Maximum number of charges on a pixel"

set label 1 "" at screen 0.03,0.45 point lt 1
set label 2 "" at screen 0.03,0.55 point lt 2
set label 3 "" at screen 0.03,0.65 point lt 3
set label 4 "" at screen 0.505,0.55 point lt 4

plot sprintf("<awk '{if($3==10.0)print}' %s | sort -k6,6n", dfname) u 6:($9/$8*2.35*100) w p t "ENC = 10 e^-/pixel",\
     sprintf("<awk '{if($3==20.0)print}' %s | sort -k6,6n", dfname) u 6:($9/$8*2.35*100) w p t "ENC = 20 e^-/pixel",\
     sprintf("<awk '{if($3==30.0)print}' %s | sort -k6,6n", dfname) u 6:($9/$8*2.35*100) w p t "ENC = 30 e^-/pixel",\
     dfname u 6:10 w p t "Maximum number of charges on a pixel" axes x1y2

unset y2label
unset y2tics
unset my2tics
set ytics mirror

set origin 0.56,-0.045
set size 0.43,1.045

set xrange [0:11]
set yrange [0:11]
set xlabel '{/Symbol s}_t [mm]'
set ylabel 'Pixel pitch [mm]'
set cbrange [0:600]
set cblabel "Number of pixels with a hit"

set palette model RGB defined (0.0 1.0 1.0 1.0, 1/256. 1.0 1.0 1.0, 1/255. 0.0 0.0 0.51, 0.34 0.0 0.81 1.0, 0.61 0.87 1.0 0.12, 0.84 1.0 0.2 0.0, 1.0 0.51 0.0 0.0) positive

set format y "%g"
plot sprintf("<awk '{if($3==10.0)print}' %s", dfname) u 1:2:6 w p pt 5 ps 2 lc palette t ''


unset multiplot

unset output
