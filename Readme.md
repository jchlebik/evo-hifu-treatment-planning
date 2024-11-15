# Evolutionary Design of Ultrasound Treatment Plans

## Introduction

This repository contains the master's thesis focusing on the evolutionary desing
of ultrasound treatment plans. The main focus is on evaluating various 
evolutionary algorithms and comparing their efficiency.


## Repository structure

    .
    +--Data       - Output data of the benchmark runs and HIFU runs, Statistical results of the same runs.
    +--Literature - Publications, references, manuals, etc.
    +--Sources    - Root folder for the sources.
    +--Thesis     - Latex sources of the thesis.
    +--Misc       - Other auxiliary materials used during the developement, such as previous work of other people.
    Readme.md     - Read me file


## Build instruction

Install basic C/C++11 build tools

Install Intel C++ Compiler

Install GAUL. Instructions and sources [here] (Sources/Thirdparty/gaul-devel-0.1850-0/)

Python script requirements are python 3.6 and SciPy

`python3 -m pip install --user numpy scipy matplotlib ipython jupyter pandas sympy nose`


## Usage instruction

To use the optimizers edit the Makefile with the appropriate fitness function, run make and provide parameters as arguments (use -h for help)

To build the GA optimizer with the default scoring function:

    `make ga`

To build the GA optimizer with a different scoring function one must edit the variable `BENCHSCOREPATH` inside the [Makefile](Sources/Makefile) so that it points to a source file containing the definitions of the desired score function in compliance with the  [FitnessFunc.h](Sources/FitnessFunc/FitnessFunc.h)
file inside [FitnessFunc folder](Sources/FitnessFunc/):

	`BENCHSCOREPATH=\....\GriewankScore.cc`
	
To build the optimizers with the HIFU system as a scoring funcion, use 

	`make ga_hifu`
	
To get help for running the GA optimizer

    `./ga -h`
  Or
  
	  `./ga_hifu -h`

PBS scripts for capacity computing in case one is running on a supporting system are provided inside the [JobScripts](Sources/JobScripts) folder.

All source code of the HIFU score funcion are inside the [HIFUScore](Sources/FitnessFunc/HIFUScore) folder. There one can also find all the data describing the medium - [AustinWoman](./FitnessFunc/HIFUScore/Data/AustinWoman/) - and also the files desciribing the target and penalization map inside [Targets](./FitnessFunc/HIFUScore/Data/Targets/). The main entry point and the definition of the entire [FitnessFunc.h](Sources/FitnessFunc/FitnessFunc.h) interface can be find inside the [heatDiffusionScore.cc](Sources/FitnessFunc/HIFUScore/heatDiffusionScore.cc) file.

To generate statistics script file [statsGenerator.py](Sources/statsGenerator.py) is provided. Use `./statsGenerator.py -h` for more detailed information on usage.

tested on Windows 10 and Ubuntu 18.04.2 LTS


## Author information

 * Name: Jakub Chlebík
 * Email: jakub.chlebik@gmail.com
 * Data: 2019/2020

