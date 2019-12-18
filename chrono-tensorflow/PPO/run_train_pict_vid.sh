#!/usr/bin/env bash

# initial training
python3 ./train_serial.py ChronoHumanoid --renderOFF -n 2000
python3 plotting_state.py 

python3 ./tester.py ChronoHumanoid --renderON --VideoSave 
# you need to rotate the figure to get a good angle for walking
echo you need to rotate the figure to get a good angle for walking
python3 combine_video.py # to compile pictures to make video 
 
echo finish


