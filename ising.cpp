#include <iostream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_integration.h>
#include <fstream>
#include <cmath>
#include "/ux/PHYSICS/PART_2/c++/cavlib/constants.hh"
#include "headers/sys_parameters.h" //include kb, mu, J, H, Ni, Nj and tmax
#include "headers/spin_ops.h"	//functions for operating on spins (flip and convert to int value)
#include "headers/calc_functions.h"	//functions to calculate various quantities of interest
#include "headers/evol_functions.h" //functions to implement a time step and initialise the spins
#include "headers/misc.h" //functions to output quantities, and a few other assorted
#include "headers/experiments.h" //functions to calculate data and print to console



int main(int argc, char* argv[]){   
  //Variables for the run
  int t_eq; //time to allow the system to reach eq. before measuring eq. quantities
  int seed_spin; //seed for selecting spins
  int seed_bolt;  //seed for boltzmann test
  int seed_init; //seed for initialization
  double Tmin;
  double Tmax;
  double dT;
  process_args(argc, argv, &t_eq, &seed_init, &seed_bolt, &seed_spin, &Tmin, &Tmax, &dT);
                                                                                                                                                  

  //Setup random number generators for spin selection and boltzmann test seperately
  const gsl_rng_type * Type;
  gsl_rng * r_init;
  gsl_rng * r_spin;
  gsl_rng * r_bolt;
  gsl_rng_env_setup();
  Type = gsl_rng_default;
  r_spin = gsl_rng_alloc (Type);
  r_bolt = gsl_rng_alloc (Type);
  r_init = gsl_rng_alloc (Type);
  gsl_rng_set(r_spin, seed_spin);
  gsl_rng_set(r_bolt, seed_bolt);
  gsl_rng_set(r_init, seed_init);



  //Declare the array of spins
  int spins[Ni][Nj];
  
  //Initialize them
  initialize_uniform(spins,1);
  //initialize_random(spins, r_init);

  //Print the parameters
  output_parameters(seed_init, seed_bolt, seed_spin, t_eq);
  
  //Calculate and print data using one of the 'experiments' 
  
  magnetisation_vs_time(spins, r_bolt, r_spin, Tmin, 'r');
  
  
  //energy_vs_time(spins, r_bolt, r_spin, Tmin, 'r');
  
  
  //snapshot(spins, t_eq, r_bolt, r_spin, Tmin);


  //magnetisation_vs_temp(Tmin, Tmax, dT, t_eq, spins, r_bolt, r_spin, 'r');
  
  
  //energy_vs_temp(Tmin, Tmax, dT, t_eq, spins, r_bolt, r_spin, 'r');
  
  
  //c_vs_temp(Tmin, Tmax, dT, t_eq, spins, r_bolt, r_spin, 'r');

  
  return 0;
}

