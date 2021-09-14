# if [ -d $1 ]
# then
   	# echo $1 " directory  exists!"

# else
   	# echo "build new directoty" $1
# fi

# mkdir $1


system=.
matter=CZT
particle=e-
entry=10000
energy1=$1

rho=5.78   # density
gdml=${matter}.gdml
input=${matter}.mac
centre=0
echo $system
echo $matter
echo $particle
echo $entry
      # ./in  ${j}
    			echo  "/control/verbose 2">${input}
			echo  "/run/verbose 2">>${input}
                        echo  "/tracking/verbose 1">>${input}
			echo  "/tracking/storeTrajectory 1">>${input}
			echo  "/G4QSim/phys/setCuts 0.001 mm">>${input}
	 		echo  "/gps/particle ${particle}">>${input}
         		#echo  "/gps/pos/type Plane">>${input}
          		#echo  "/gps/pos/shape Circle">>${input}
          		#echo  "/gps/pos/rot1 1.0 0. 0.">>${input}
          		#echo  "/gps/pos/centre 0. 0. 0.">>${input}
          		#echo  "/gps/pos/radius 1.5 mm">>${input}
          		#echo  "/gps/ang/type iso">>${input}
			echo  "/gps/energy ${energy1} keV">>${input}
			echo  "/gps/position 0.0 0.0 0.0 mm">>${input}
			echo  "/gps/direction 0.0 0.0 1.0">>${input}


	./../../data_generator/geant4.9.6.2/Linux-g++/G4QSim -g ${gdml} -o ${system}/e.root -n ${entry} -m ${input} >${system}/log_e.txt
	rm ${system}/log_e.txt
#	root -l -q -b "${system}/nengsunlv.C("\"CZT_energy_deposited.txt\"")"

