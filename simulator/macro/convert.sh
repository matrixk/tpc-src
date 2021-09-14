#!/bin/bash

dataFile=$1
lineCount=$(awk '{print NR}' $dataFile|tail -n1)
eventNumber=`expr $lineCount / 3`
echo $dataFile $lineCount $eventNumber

touch new.dat

for i in $( seq 1 $eventNumber)
do
	NUM=$i
	l1=`expr 3 \* $NUM - 2`
	l2=`expr 3 \* $NUM - 1`
	l3=`expr 3 \* $NUM`
	type1=$(awk 'NR=='$l1' {printf  $1}' $dataFile)
	echo "       "$NUM"  .000000       "$type1 >> new.dat
	type2_1=$(awk 'NR=='$l2' {printf  $1}' $dataFile)
	p1_1=$(awk 'NR=='$l2' {printf "%.6f\n\n",  $2}' $dataFile)
	p2_1=$(awk 'NR=='$l2' {printf "%.6f\n",  $3}' $dataFile)
	p3_1=$(awk 'NR=='$l2' {printf "%.6f\n",  $4}' $dataFile)
	p4_1=$(awk 'NR=='$l2' {printf "%.6f\n",  $5}' $dataFile)
	echo "  "$type2_1"  "$p2_1"      "$p3_1"      "$p4_1"      "$p1_1  >> new.dat
	type2_2=$(awk 'NR=='$l3' {printf  $1}' $dataFile)
        p1_2=$(awk 'NR=='$l3' {printf "%.6f\n",  $2}' $dataFile)
        p2_2=$(awk 'NR=='$l3' {printf "%.6f\n",  $3}' $dataFile)
        p3_2=$(awk 'NR=='$l3' {printf "%.6f\n",  $4}' $dataFile)
        p4_2=$(awk 'NR=='$l3' {printf "%.6f\n",  $5}' $dataFile)
        echo "  "$type2_2"  "$p2_2"      "$p3_2"      "$p4_2"      "$p1_2  >> new.dat
	echo "Done event "$NUM
done

