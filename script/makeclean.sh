# if [ -d $1 ]
# then
   	# echo $1 " directory  exists!"
# else
   	# echo "build new directoty" $1
# fi
# mkdir $1
chmod 777 *.sh
cd ../data_generator
make clean

cd ../analysis/energy_loss
make clean
cd ../diffusion
make clean
cd ../sigma_over_sqrt_de
make clean

cd ../../script
