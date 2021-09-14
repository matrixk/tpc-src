#!/bin/bash


# detector matter 
matter=CZT

# particle=neutron
#particle list
#particle=e+
#particle=e-
particle=gamma
#particle=proton
#particle=neutron

#how many events
entry=100

# particle energy in KeV and Z axis, energy can scan
energy1=500
energy2=500
energy_gap=50

# matter thickness in Z and X Y, 
thick=100
# thick can scan in list
thickx=400
thicky=400
# thick is in mm
centre_control=1
# the material centre is 0 0 thick/2 if centre_control=1
# the material centre is 0 0 0       if centre_control=0

rho=1   # density not used; only appear in data file name


#parameters for analysis plot ^
thick_start=0
thick_stop=100
thick_step=1
# thick is in mm
diffusion_range=10
diffusion_xybin=600
diffusion_zbin=150
digitization_output_name=digitization_output_name.txt
#parameters for analysis plot v


#################################################
system=../data
gdml=${system}/${matter}.gdml
input=${system}/${matter}.mac
centre=0
echo $system
echo $matter
echo $particle
echo $entry
for((i=${energy1};i<=${energy2};i+=${energy_gap}))
#  for i in 20 40 60 80 100 200 300 400 500 600 700 800
   do
       echo $i
  for j in $rho
      do

     for houdu in ${thick}
          do
      # ./in  ${j}
    			echo  "/control/verbose 0">${input}
			echo  "/run/verbose 0">>${input}
			echo  "/G4QSim/phys/addPhysics penelope">>${input}
			echo  "/G4QSim/phys/addPhysics neutronHP">>${input}
			# echo  "/G4QSim/phys/addPhysics gammanuclear">>${input}
			# echo  "/G4QSim/phys/addPhysics ionphysics">>${input}
                        echo  "/tracking/verbose 0">>${input}
			echo  "/tracking/storeTrajectory 1">>${input}
			echo "/run/initialize">>${input}
			echo  "/G4QSim/phys/setCuts 0.001 mm">>${input}
			#echo  "/G4QSim/phys/setGCut 0.001 mm">>${input}
			#echo  "/G4QSim/phys/setECut 0.001 mm">>${input}
			#echo  "/G4QSim/phys/setPCut 0.001 mm">>${input}
			#echo  "/G4QSim/phys/setProtonCut 0.1 mm">>${input}
			#echo  "/G4QSim/phys/addPhysics neutronHP">>${input}
			#echo  "/testhadr/phys/ QGSP">>${input}
	 		echo  "/gps/particle ${particle}">>${input}
         		#echo  "/gps/pos/type Plane">>${input}
          		#echo  "/gps/pos/shape Circle">>${input}
          		#echo  "/gps/pos/rot1 1.0 0. 0.">>${input}
          		#echo  "/gps/pos/centre 0. 0. 0.">>${input}
          		#echo  "/gps/pos/radius 1.5 mm">>${input}
          		#echo  "/gps/ang/type iso">>${input}
			echo  "/gps/energy ${i}.0 keV">>${input}
			echo  "/gps/position 0.0 0.0 0.0 mm">>${input}
			echo  "/gps/direction 0.0 0.0 1.0">>${input}
#        ./inc CZT_initial.gdml CZT.gdml 64 ${houdu}
#        ./inc initial.gdml ${matter}.gdml 64 ${houdu}
	echo ${houdu}
	
if [ ${centre_control} != 1 ]
then
   	centre=0
else
   	centre=$(echo "${houdu}/2"|bc -l)
fi
    
	data_file=${system}/${matter}_${particle}_E${i}_thick${houdu}_rho${j}.root
		
		echo centre= ${centre}
		cp ./../geometry/1.gdml  ${system}/1.gdml
		echo   "  <position name=\"gasp\" unit=\"mm\" x=\"0\" y=\"0\" z=\"${centre}\" />" >> ${system}/1.gdml
        cat ${system}/1.gdml ./../geometry/2.gdml > ${system}/2.gdml 
        echo   "  <box aunit=\"radian\" lunit=\"mm\" name=\"gasV\" x=\"${thickx}\" y=\"${thicky}\" z=\"${houdu}\" />" >> ${system}/2.gdml
		cat ${system}/2.gdml ./../geometry/3.gdml > ${system}/3.gdml
        echo   "   <materialref ref=\"${matter}\" />" >> ${system}/3.gdml
		cat ${system}/3.gdml ./../geometry/4.gdml > ${system}/${matter}.gdml


  
if [ x$1 != xana ]
then
   	./../data_generator/Linux-g++/G4QSim -g ${gdml} -o ${data_file} -n ${entry} -m ${input} 
else
   	echo analysis only
fi


#################### analysis

if [ x$1 != xdata ]
then

echo  ${data_file} >data_file.txt
echo  ${digitization_output_name} >data_out.txt

# ./../analysis/energy_loss/energy_loss -input ${data_file} -output energy_loss.root -start ${thick_start} -stop ${thick_stop} -step ${thick_step} -energy ${i} -energy_bin ${i} -gr_name de_z -ratio_name de/e_z


# ./../analysis/diffusion/diffusion -input ${data_file} -output diffusion.root -start ${thick_start} -stop ${thick_stop} -step ${thick_step} -diffusion_range ${diffusion_range} -diffusion_xybin ${diffusion_xybin} -diffusion_zbin ${diffusion_zbin} 

# ./../analysis/sigma_over_sqrt_de/sigma_over_sqrt_de -input_de energy_loss.root -input_sigma diffusion.root -output sigma_over_sqrt_de.root -start ${thick_start} -stop ${thick_stop} -step ${thick_step} 

# ./../analysis/digitization/digitization -input ${data_file} -output_name ${digitization_output_name} -Amplification 100 -Threshold 50
./../analysis/pixelSizeTradeOff/geant2root  para.txt

# ./../analysis/cluster/cluster -input_name ${digitization_output_name} -photo_on $2

	
else
   	echo data generation only
fi


	done
	done
        done
