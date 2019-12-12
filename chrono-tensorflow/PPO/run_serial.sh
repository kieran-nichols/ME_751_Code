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
python3 ./train_serial.py ChronoAnt --renderOFF -b 1 -n 1
