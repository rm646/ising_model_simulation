#ifndef EVOL_FUNCTIONS_H
#define EVOL_FUNCTIONS_H

//Function declaration and description
//====================================
//functions to evolve the system
void step_rand(int spins[Ni][Nj], gsl_rng* r_bolt, gsl_rng* r_spin, const double T); //execute a time step by randomly selecting Ni*Nj spins using the rng specified
void step_syst(int spins[Ni][Nj], gsl_rng* r_bolt, const double T); //selects all Ni*Nj spins systematically
bool boltzmann(double dE, double* p, const double T); //tests whether exp(-dE/kT)>p where p is [0,1] and drawn from a uniform dist.
void initialize_uniform(int spins[Ni][Nj], int up); //uniform initialisation of spins, up_or_down=1 if up desired
void initialize_random(int spins[Ni][Nj], gsl_rng* r); //random initialisation


//Function implementations
//====================================
void step_rand(int spins[Ni][Nj], gsl_rng* r_bolt, gsl_rng* r_spin, const double T){
  int max = Ni*Nj;
  int i=0;
  int j=0;
  double dE=-64; //initialise to unusual value for error checking
  double p; //to hold the random number used to compare the Boltzmann factor to
  for(int k=0; k<max; k++){
    i = gsl_rng_uniform_int(r_spin, Ni);
    j = gsl_rng_uniform_int(r_spin, Nj);
    
    dE = calc_dE(spins, i, j);

    p = gsl_rng_uniform(r_bolt);
    if (dE<=0)
      {flip(spins[i][j]);}
      
    else if (boltzmann(dE, &p, T))
      {flip(spins[i][j]);}    
  }
}



void step_syst(int spins[Ni][Nj], gsl_rng* r_bolt, const double T){
  //Steps through the array of spins, testing for flipping using boltzmann() function and <0 condition
  double dE;
  double p;
  for(int i=0; i<Ni; i++){
    for(int j=0; j<Nj; j++){
    
      dE = calc_dE(spins, i, j);
      
      p = gsl_rng_uniform(r_bolt);
      if(dE<=0)
	{flip(spins[i][j]);}
      else{ if (boltzmann(dE, &p, T)){flip(spins[i][j]);}}
    }
  }
}



bool boltzmann(double dE, double* p, const double T){
  if(gsl_sf_exp(-dE/T)>*p){return true;}
  else{return false;}
}




void initialize_uniform(int spins[Ni][Nj], int spin){
//if spin is not an acceptable value, spins just remain uninitialized
  if (spin==1){
    for(int i=0; i<Ni; i++){
      for(int j=0; j<Nj; j++){
	spins[i][j]=spin;
      }
    }
  }
  if (spin==-1){
    for(int i=0; i<Ni; i++){
      for(int j=0; j<Nj; j++){
	spins[i][j]=spin;
      }
    }
  }
}



void initialize_random(int spins[Ni][Nj], gsl_rng* r){
  for(int i=0; i<Ni; i++){
    for(int j=0; j<Nj; j++){
      spins[i][j] = 2*gsl_rng_uniform_int(r,2) -1;
    }
  }
}



void output_spins(int spins[Ni][Nj]){
  for(int i=0; i<Ni; i++){
    for(int j=0; j<Nj; j++){
      std::cout<<i<<"\t"<<j<<"\t"<<spins[i][j]<<"\n";
    }
    std::cout<<"\n";
  }
}



#endif
