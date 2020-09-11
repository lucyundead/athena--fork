#!/bin/bash

make clean

# for isothermal
#python ./configure.py --prob bar --coord cartesian --eos isothermal --nghost 3 --flux roe --cxx icc -mpi -hdf5 --hdf5_path /share/home/jtshen/soft/hdf5


python ./configure.py --prob milkyway --coord cartesian --eos isothermal --nghost 2 --flux roe --cxx icc -mpi -hdf5 --hdf5_path /share/home/jtshen/soft/hdf5

# for adiabatic + cooling
#python ./configure.py --prob bar --coord cartesian --eos adiabatic  --nghost 3 --flux hllc --cxx icc -mpi -hdf5 --hdf5_path /share/home/jtshen/soft/hdf5


#python ./configure.py --prob milkyway --coord cartesian --eos adiabatic --nghost 2 --flux hllc --cxx icc -mpi -hdf5 --hdf5_path /share/home/jtshen/soft/hdf5


make

