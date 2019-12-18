# Machine Learning using PyChrono and Tensorflow

## Getting Started

### Prerequisites
PyChrono/Chrono,
Tensorflow,
Bash,
Python

### Installing
Bash: Use linux system (this project was tested on Ubuntu)

Python: Anaconda (conda)

PyChrono: Refer to http://projectchrono.org/pychrono/

Tensorflow:Refer to https://docs.anaconda.com/anaconda/user-guide/tasks/tensorflow/


## Project Description
My project centers about utilizing PyChrono and Machine Learning to create a controller that achieves walking for the multi-segmented body of a four-legged ant and humanoid.  

### Goals
Implementation of Gym-Chrono: machine learning component of PyChrono on personal computer with visualization working  
Scaling analysis to show differences of personal computer using serial vs parallel code  
Creation of humanoid environment and report of success/failure with alterations and additions of various cost functions  

### Environments
ChronoAnt
Picture
ChronoHumanoid
Picture

### Machine Learning Algorithm
The machine learning for this project uses coding structures from Tensorflow. Tensorflow is an open-source machine learning platform that allows for easy model building and production of robust machine learning. To solve tasks, neural networks will be used in the process of Proximal Policy Optimization (PPO), which is a form of reinforcement learning. Reinforcement learning is the process promoting actions of a program in order to maximize its specified cumulative reward. PPO further specifies that as the machine learning be based in a stochastic gradient ascent (way of maximizing reward) to perform each set of iterations (in my case; ensuring that the robot stays up).


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
