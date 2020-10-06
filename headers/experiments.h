#ifndef EXPERIMENTS_H
#define EXPERIMENTS_H

//Function declaration and description
//====================================
//Calculates and outputs (to console) the spatial average of M (per site) for each time step. option may be either s or r for 
//systematic or random selection of spins for probing respectively.
//Outputs: t "\t" M
void magnetisation_vs_time(int spins[Ni][Nj], gsl_rng* r_bolt, gsl_rng* r_spin, const double T, char option);


//Calculates and outputs (to console) the spatial average of M (per site) for each time step. option may be either s or r for 
//systematic or random selection of spins for probing respectively.
//Outputs: t "\t" E
void energy_vs_time(int spins[Ni][Nj], gsl_rng* r_bolt, gsl_rng* r_spin, const double T, char option);


//Calculates and outputs to console the time average (over the period (tmax-t_eq)) of magnetisation
//And the standard dev. for a series of temperatures
//Arguments: Tmin, Tmax, dT, t_eq, spins, r_bolt, r_spin, option.
//Outputs: T "\t" <M> "\t" devM
void magnetisation_vs_temp(double Tmin, double Tmax, double dT, int t_eq, int spins[Ni][Nj], gsl_rng* r_bolt, gsl_rng* r_spin, char option);


//Similar to function for magnetisation vs temp
//Outputs: T "\t" <E> "\t" devE
void energy_vs_temp(double Tmin, double Tmax, double dT, int t_eq, int spins[Ni][Nj], gsl_rng* r_bolt, gsl_rng* r_spin, char option);


//Calculates and outputs the heat capacity as estimated from dE/dT and fluctuation dissipation theorem
//Arguments: Tmin, Tmax, dT, t_eq, spins, r_bolt, r_spin, option
//Outputs T, C(discrete method), C(fluctuation-diss), C(analytic), error in integral
void c_vs_temp(double Tmin, double Tmax, double dT, int t_eq, int spins[Ni][Nj], gsl_rng* r_bolt, gsl_rng* r_spin, char option);



//Function implementations
//================================
void magnetisation_vs_time(int spins[Ni][Nj], gsl_rng* r_bolt, gsl_rng* r_spin, const double T, char option){
  if(option!='r' && option !='s'){std::cerr<<"option not valid\n\n"; return;}
  std::cout<<"#Experiment: magnetisation vs time\n"<<"#t\tmag_per_spin\n";
  if(option=='s'){
    double N = (double)(Ni*Nj);
    for(int t=0; t<tmax; t++){
      step_syst(spins, r_bolt, T);
      std::cout<<t<<"\t"<<( mu*((double)calc_total_S(spins)) )/N<<"\n";
    }
  }
	
	
  if(option=='r'){
    double N = (double)(Ni*Nj);
    for(int t=0; t<tmax; t++){
      step_rand(spins, r_bolt, r_spin, T);
      std::cout<<t<<"\t"<<( mu*((double)calc_total_S(spins)) )/N<<"\n";
    }
  }
}



void energy_vs_time(int spins[Ni][Nj], gsl_rng* r_bolt, gsl_rng* r_spin, const double T, char option){
  if(option!='r' && option !='s'){std::cerr<<"option not valid\n\n"; return;}
  std::cout<<"#Experiment: energy vs time\n"<<"#t\tE_per_spin\n";
  if(option=='s'){
    double N = (double)(Ni*Nj);
    for(int t=0; t<tmax; t++){
      step_syst(spins, r_bolt, T);
      std::cout<<t<<"\t"<<(calc_total_E(spins)/N)<<"\n";
    }
  }
	
	
  if(option=='r'){
    double N = (double)(Ni*Nj);
    for(int t=0; t<tmax; t++){
      step_rand(spins, r_bolt, r_spin, T);
      std::cout<<t<<"\t"<<(calc_total_E(spins)/N)<<"\n";
    }
  }
}


void magnetisation_vs_temp(double Tmin, double Tmax, double dT, int t_eq, int spins[Ni][Nj], gsl_rng* r_bolt, gsl_rng* r_spin, char option){
  if(option!='r' && option!='s'){std::cerr<<"option not valid\n\n"; return;}
  std::cout<<"#Experiment: magnetisation vs temp\n"<<"#T\tav_mag_per_spin\terror_in_it\tanalytic_value\n";
	
  //Set up variables to hold results		
  double N = (double)(Ni*Nj);
  int S[tmax]; //array of total spin at various times
  double avS;
  double devS;
  double analytic=0;
  		
  int n=0;
  double dS = 0; //lets us seek the largest change, to then find the heat capacity
  double dSmaxsq=0;
  double Tcrit=0;
  		
  //double anaerr=0;
  	
  //Vary the temperature as specified by the arguments	
  for(double T=Tmin; T<Tmax; T+=dT){
		
    if(option=='s'){
      for(int t=0; t<tmax; t++){
	step_syst(spins, r_bolt, T);
	if(t>=t_eq){S[t] = calc_total_S(spins);}				
      }
    }
			
    if(option=='r'){
      for(int t=0; t<tmax; t++){
	step_rand(spins, r_bolt, r_spin, T);
	if(t>=t_eq){S[t] = calc_total_S(spins);}				
      }
    }	
				
    if(n>0){dS=avS;}
			
    calc_av_S(S, avS, devS, t_eq);
			
    dS -= avS;
    dS = dS*dS;
		       
    if(n>0){if(dS>dSmaxsq){dSmaxsq = dS; Tcrit = T-(0.5*dT);}}
			
    analytic = calc_exact_M_perspin(T);
    std::cout<<T<<"\t"<<mu*(avS/N)<<"\t"<<mu*(devS/N)<<"\t"<<analytic<<"\t"<<dS<<"\n";
    n++;
  }
  std::cout<<"\n"<<Tcrit<<"\n";	
	
}



void energy_vs_temp(double Tmin, double Tmax, double dT, int t_eq, int spins[Ni][Nj], gsl_rng* r_bolt, gsl_rng* r_spin, char option){
  if(option!='r' && option!='s'){std::cerr<<"option not valid\n\n"; return;}
  std::cout<<"#Experiment: energy vs temp\n"<<"#T\tav_E_per_spin\terror_in_it\tanalytic_value\ttheory_error\n";	
	
  //Set up variables to hold results		
  double N = (double)(Ni*Nj);
  double E[tmax]; //array of total spin at various times
  double avE;
  double devE;
  double analytic=0;
  double anaerr=0;
  
  int n=0;
  double dE = 0;
  double dEmaxsq =0;
  double Tcrit=0;
  		
		

  //Vary the temperature as specified by the arguments	
  for(double T=Tmin; T<Tmax; T+=dT){
    if(option=='s'){
      for(int t=0; t<tmax; t++){
	step_syst(spins, r_bolt, T);
	if(t>=t_eq){E[t] = calc_total_E(spins);}				
					
      }
    }
			
    if(option=='r'){
      for(int t=0; t<tmax; t++){
	step_rand(spins, r_bolt, r_spin, T);
	if(t>=t_eq){E[t] = calc_total_E(spins);}				
      }
    }
			
    if(n>0){dE=avE;}
		
    calc_av_E(E, avE, devE, t_eq);
			
    dE -= avE;
    dE  = dE*dE;
    if(n>0){if(dE>dEmaxsq){dEmaxsq = dE; Tcrit = T-(0.5*dT);}}
			
    analytic = calc_exact_E_perspin(T, &anaerr);
    std::cout<<T<<"\t"<<(avE/N)<<"\t"<<(devE/N)<<"\t"<<analytic<<"\t"<<anaerr<<"\t"<<dE<<"\n";
    n++;
  }	
  std::cout<<"\n"<<Tcrit<<"\n";	
}



void c_vs_temp(double Tmin, double Tmax, double dT, int t_eq, int spins[Ni][Nj], gsl_rng* r_bolt, gsl_rng* r_spin, char option){
  if(option!='r' && option!='s'){std::cerr<<"option not valid\n\n"; return;}
  std::cout<<"#Experiment: C vs temp (discrete dE/dT method)\n"<<"#(T-dT)<<\t<<((dE/dT)/N)<<\t<<Cfluc<<\t<<analytic<<\t<<anaerr<<\t<<anaest\n";		
  //Set up variables to hold results		
  double E[tmax]; //array of total spin at various times
  double avE;
  double devE;
  double dE;
  double Cfluc;
  double Cdisc;
  double analytic=0;
  double anaerr=0;
  
  int n=0;
  double N = (double)(Ni*Nj);
		
  double Cmax_disc=0;	//for estimates of Tcrit for each method
  double Cmax_fluc=0;
  double Tcrit_disc=0;
  double Tcrit_fluc=0;
		
  //Vary the temperature as specified by the arguments	
  for(double T=Tmin; T<Tmax; T+=dT){
    if(option=='s'){ 
      for(int t=0; t<tmax; t++){
	step_syst(spins, r_bolt, T);
	if(t>=t_eq){E[t] = calc_total_E(spins);}				
      }		
    }
    if(option=='r'){
      for(int t=0; t<tmax; t++){
	step_rand(spins, r_bolt, r_spin, T);
	if(t>=t_eq){E[t] = calc_total_E(spins);}				
      }	
    }
		   
    if(n>0){
      Cfluc = (devE*devE)/(Ni*Nj*(T-dT)*(T-dT)); 
      if(Cfluc*Cfluc>Cmax_fluc){Cmax_fluc=Cfluc; Tcrit_fluc= T-dT;}
    }
		  
    dE = avE;
    calc_av_E(E, avE, devE, t_eq);
    dE = avE - dE;
    Cdisc = ((dE/dT)/N);
    if(n>0){
      if(Cdisc*Cdisc>Cmax_disc){Cmax_disc=Cdisc; Tcrit_disc= T-dT;}
    }
		  
    analytic = calc_exact_c_perspin((T-dT), &anaerr);
    
    if(n>0){std::cout<<(T-dT)<<"\t"<<Cdisc<<"\t"<<Cfluc<<"\t"<<analytic<<"\t"<<anaerr<<"\n";}
    n++;
  }
  std::cout<<"\n#Fluc:\t"<<Tcrit_fluc<<"\tDisc:\t"<<Tcrit_disc<<"\n";	
}


#endif
