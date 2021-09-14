# if [ -d $1 ]
# then
   	# echo $1 " directory  exists!"
# else
   	# echo "build new directoty" $1
# fi
# mkdir $1

#chmod 777 *.sh
#cd ../data_generator
#make

#cd ../analysis/energy_loss
#make
#cd ../diffusion
#make
#cd ../sigma_over_sqrt_de
#make



chmod 777 *.sh
cd ../data_generator
make

cd ../analysis/energy_loss
rm energy_loss
make


cd ../diffusion
rm diffusion
make

cd ../sigma_over_sqrt_de
rm sigma_over_sqrt_de
make

cd ../file_print
rm file_print
make

cd ../digitization
rm digitization
make

cd ../cluster
rm cluster
make


cd ../read_data
rm read_data
make

cd ../../script
