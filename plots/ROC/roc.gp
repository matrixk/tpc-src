dfname="roc.dat"
dsets="10x10x4 2x2x2 2x2x4 4x4x2 4x4x4 6x6x2 6x6x4 8x8x2 8x8x4 3D:8x8x8"
dnames="xy10z4 xy2z2 xy2z4 xy4z2 xy4z4 xy6z2 xy6z4 xy8z2 xy8z4 3D:xy8z8"
dlist="1 2 3 4 5 6 7 8 0 9"

set output "1.eps"
set terminal postscript eps enhanced solid color "Helvetica" 20
#set terminal postscript eps enhanced dl 2.0 lw 2.0 size 3.25,2.275 "Helvetica" 14

set multiplot

set size 0.78, 0.75
set origin 0.1, 0.15

set key left bottom

unset xlabel
set xrange [0:0.85]
unset ylabel
set yrange [0.9:1.001]
set ytics 0.01

set y2tics ("10" 0.9, "20" 0.95, "100" 0.99)
set y2label "Suppression factor"

set grid

plot for [i in dlist] dfname index (i+0) u (1.0-$1):2 w l lw ((i+1)*0.8) t word(dnames, (i+1))

set size 1,1
set origin 0,0

unset y2tics
unset y2label

unset log y
set ytics autofreq
set ylabel 'Background rejection rate'
set yrange [0.0:1.01]
set xlabel 'Signal efficiency'
set xrange [0.0:1.0]

unset grid

plot for [i in dlist] dfname index (i+0) u (1.0-$1):2 w l t ''

unset multiplot

unset output
set term x11
