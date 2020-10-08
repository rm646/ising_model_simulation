
## Introduction

I wrote this code sometime in 2013 as part of my undergraduate course.
I'm sharing it now as demonstration of my C++ and modelling experience,
and also to embarass me into better code structure.

## Setup

[Gnu Scientific Library (GSL)](https://www.gnu.org/software/gsl/) 1.16 is
required to run this.

I compiled it with
```
g++ -O3 -Wall ising.cpp -o ising -lgsl
```

### Usage

'ising.cpp' performs the initialization from command line arguments. Different
'experiments' (from headers/experiments.h) are chosen by commenting out the
others in 'ising.cpp'. The program will output line by line to the console.

```
./ising <t_eq> <seed_init> <seed_bolt> <seed_spin> <Tmin> <Tmax> <dT>
```

Arguments:
- t_eq: time to allow the system to reach eq. before measuring eq. quantities
- seed_spin: seed for selecting spins
- seed_bolt:  seed for boltzmann test
- seed_init: seed for initialization
- Tmin: the set temperature for experiments which don't vary in temperature, and the minimum temperature for experiments which do
- Tmax: the maximum temperature for experiments which vary in temperature
- dT: the temperature step for experiments which vary in temperature
