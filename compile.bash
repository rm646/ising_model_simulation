#!/bin/bash
source /ux/PHYSICS/PART_2/c++/bashrc
g++ -O3 -pg -Wall $PT2_INC ising.cpp -o ising $PT2_LIB
