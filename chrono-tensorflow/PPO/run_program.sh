#!/usr/bin/env bash

###################################
## commented code for use on euler; didn't get to work
#SBATCH --job-name=test
#SBATCH -c 1
#SBATCH --time=0-0:01:00
#SBATCH --output="Output.txt"

#module load anaconda/3
#bootstrap_conda
#conda create --name virtualenv python=3.7
#conda activate virtualenv
#module load irrlicht/1.8.3
#conda install numpy
#conda install tensorflow
#conda install scikit-learn
###################################

# test
#python3 ./tester.py ChronoAnt --renderOFF -b 1 -n 1

# scaling anlysis
end=3
>Scaling.csv
for ((i=1;i<end;i++))
	do
	num=$((10 **($i)))
	echo -n "$num ">>Scaling.csv
	printf ", ">> Scaling.csv
	
	start=$seconds
	python3 ./train_serial.py ChronoAnt --renderOFF -n $num
	duration=$(( SECONDS - start ))
	echo -n "$duration ">>Scaling.csv;
	printf ", ">> Scaling.csv
	
	# plot state at last iteration of loop for the serial code
	if ((i==$end-1))
		then
		python3 plotting_state.py 
	fi
	
	start=$seconds
	python3 ./train_parallel.py ChronoAnt -n $num
	duration=$(( SECONDS - start ))
	echo -n "$duration ">>Scaling.csv
	
	printf "\n">> Scaling.csv
	echo progress: $i out of $(( end-1 ))
	done 
python3 plotting_scaling.py 
echo finished scaling with variable total iteration number

# Scaling with variable batch sizes
>Scaling_b.csv
itr_size=100
for ((i=1;i<end;i++))
	do
	num=$((10 **($i)))
	echo -n "$num ">>Scaling_b.csv
	printf ", ">> Scaling_b.csv
	
	start=$seconds
	python3 ./train_serial.py ChronoAnt --renderOFF -n $itr_size -b $num
	duration=$(( SECONDS - start ))
	echo -n "$duration ">>Scaling_b.csv;
	printf ", ">> Scaling_b.csv
	
	# plot state at last iteration of loop for the serial code
	#if ((i==$end-1))
	#	then
	#	python3 plotting_state.py 
	#fi
	
	start=$seconds
	python3 ./train_parallel.py ChronoAnt -n $itr_size -b $num
	duration=$(( SECONDS - start ))
	echo -n "$duration ">>Scaling_b.csv
	
	printf "\n">> Scaling_b.csv
	echo progress: $i out of $(( end-1 ))
	done 
python3 plotting_scaling_b.py 
echo finished scaling with variable batch numbers

