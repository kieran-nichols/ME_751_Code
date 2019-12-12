#!/usr/bin/env bash
#SBATCH --job-name=test
#SBATCH -c 1
#SBATCH --time=0-0:01:00
#SBATCH --output="Output.txt"

module load anaconda/3
bootstrap_conda
#conda create --name virtualenv python=3.7
conda activate virtualenv
module load irrlicht/1.8.3
#conda install numpy
#conda install tensorflow
#conda install scikit-learn

#python3 ./tester.py ChronoAnt --renderOFF -b 1 -n 1

end=10
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
	
	start=$seconds
	python3 ./train_parallel.py ChronoAnt -n $num
	duration=$(( SECONDS - start ))
	echo -n "$duration ">>Scaling.csv
	
	printf "\n">> Scaling.csv
	echo progress: $i out of $(( end-1 ))
	done 
echo finished

python3 plotting_scaling.py 
