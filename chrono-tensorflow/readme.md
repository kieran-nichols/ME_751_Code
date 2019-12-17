# Machine Learning using PyChrono and Tensorflow

## Getting Started

### Prerequisites
PyChrono/Chrono
Tensorflow
Bash
Python

### Installing
Bash
Use linux system (tested on Ubuntu

Python
Anaconda (conda)

PyChrono
Refer to http://projectchrono.org/pychrono/

Tensorflow
Refer to https://docs.anaconda.com/anaconda/user-guide/tasks/tensorflow/

Gym chrono?

## Project Description
My projects centers about utilizing PyChrono and Machine Learning to create a controller that achieves walking for the multi-segmented body of a four-legged ant and humanoid.  

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

`python ./train_serial.py ChronoAnt -n 20000 -b 100  --renderOFF`
|---|---|---|---|---|---|
| 1	     | 2	            | 3	       | 4      | 5     |  6       |
|---|---|---|---|---|---|

1. Use python or python3 command in a bash terminal
2. This command specifies that the serial learning will be utilized. './train_parall.py' can also be used for parallel learning
3. 'ChronoAnt' will use the Ant object. 'ChronoHumanoid' will use the humanoid object
4.
5.
6.

Other commands
Discount factor
Lambda
Kullback Leibler divergence target value


### Output

Checkpoints for saving

Video saving

Video Playing

Sacling analysis

### Problems
Couldn't get it to work on Windows or WSL on Windows  
Gymchrono vs Pychrono  
Ending execution  

## Other

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
