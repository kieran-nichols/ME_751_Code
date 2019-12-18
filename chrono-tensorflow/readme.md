# Machine Learning using PyChrono and Tensorflow

## Getting Started

### Prerequisites
PyChrono/Chrono
Tensorflow
Bash
Python

### Installing
Bash
Use linux system (this project was tested on Ubuntu)

Python
Anaconda (conda)

PyChrono
Refer to http://projectchrono.org/pychrono/

Tensorflow
Refer to https://docs.anaconda.com/anaconda/user-guide/tasks/tensorflow/


## Project Description
My project centers about utilizing PyChrono and Machine Learning to create a controller that achieves walking for the multi-segmented body of a four-legged ant and humanoid.  

### Goals
Implementation of Gym-Chrono: the machine learning component of PyChrono on personal computer with visualization working  
Scaling analysis to show differences of personal computer using serial vs parallel code  
Creation of humanoid environment and report of success/failure with alterations and additions of various cost functions  

### Environments
ChronoAnt
Picture
ChronoHumanoid
Picture

### Machine Learning Algorithm
PPO


### Program Execution

`python ./train_serial.py ChronoAnt -n 20000 -b 100  --renderON`

1. Use python or python3 command in a bash terminal
2. This command specifies that the serial learning will be utilized. `./train_parallel.py` can also be used for parallel learning. `./tester.py` can be used for running the episodes without further improving the learning. 
3. 'ChronoAnt' will use the Ant robot. 'ChronoHumanoid' will use the humanoid robot. The robot are also called env (environment)
4. `-n 20000` says that 20000 episodes will be used
5. `-b 100` says that batches of 100 will be processed at a time for learning.
6. `--renderON` displays the robot and floor as it learns how to walk or fails to walk.

Other commands
Discount factor: '-g'
Lambda: '-l'
Kullback Leibler divergence target value: '-k'
Help: '-h'


### Output

Checkpoints for saving

Video saving

Video Playing

Sacling analysis

### Problems
Couldn't get it to work on Windows or WSL on Windows  
There are various versions of Machine learning for PyChrono available and it's difficult to choose.


## Author

[Kieran Nichols](https://github.com/kieran-nichols)

## License

## Acknowledgements

[Chrono](http://projectchrono.org/),  [PyChrono](https://projectchrono.org/pychrono/),  [Simone Benatti](https://github.com/Benatti1991),  [Dan Negrut](http://homepages.cae.wisc.edu/~negrut/),  [Keanu Nichols](https://github.com/kmn5409)

## References
https://www.markdownguide.org/cheat-sheet/
https://github.com/projectchrono/gym-chrono
https://github.com/projectchrono/chrono/tree/develop/src/demos/python/chrono-tensorflow
http://api.projectchrono.org/tutorial_pychrono_demo_tensorflow.html
