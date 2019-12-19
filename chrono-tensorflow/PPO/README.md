# PPO

## Program Execution
Basic program execution

`python ./train_serial.py ChronoAnt -n 20000 -b 100  --renderON`

1. Use python or python3 command in a bash terminal
2. This command specifies that the serial learning will be utilized. `./train_parallel.py` can also be used for parallel learning. `./tester.py` can be used for running the episodes without further improving the learning. 
3. 'ChronoAnt' will use the Ant robot. 'ChronoHumanoid' will use the humanoid robot. The robot are also called env (environment)
4. `-n 20000` says that 20000 episodes will be used
5. `-b 100` says that batches of 100 will be processed at a time for learning.
6. `--renderON` displays the robot and floor as it learns how to walk or fails to walk.

Other commands

Discount factor: `-g`,
Lambda: `-l`,
Kullback Leibler divergence target value: `-k`,
Help: `-h`

## Project's high level program
This project has a shell script that can be run in bash to have outputs that log checkpoints for saving, save pictures of rendering, processing videos from pictures, scaling analysis of the serial and parallel codes, and create state graphs that describe the position and velocity of the robots.

To execute program, use `bash run_program.sh`

`Ctrl-C` doesn't usually work for stopping the python/bash program when using PyChrono and TensorFlow. `pkill -9 python` will kill the python kernel.
 
After running an individual commands of 
`python3 ./train_serial.py ChronoHumanoid --renderOFF -n 2000`
`python3 plotting_state.py `
`python3 ./tester.py ChronoHumanoid --renderON –VideoSave`
`python3 combine_video.py # to compile pictures to make video`

You can run `run_train_pict_vid.sh` but you should read the code before running. You will need to rotate the rendering to get a good camera view of the robot.