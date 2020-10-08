#ifndef MISC_H
#define MISC_H

//Function declaration and description
//====================================
//functions to ouput various things to the console
void output_spins(int spins[Ni][Nj]); //outputs x index, y index, and spin value
void output_spin_img(int spins[Ni][Nj]); //ouputs the spins as a grid of ASCII art
void output_parameters(int seed_init, int seed_bolt, int seed_spin, int t_eq); //outputs all parameters to console, prefixed by #


//function to sort and set parameters using command line arguments
int process_args(int argc, char* argv[], int* t_eq, int* seed_init, int* seed_bolt, int* seed_spin, double* T);

//Function implementations
//====================================

//function to take care of command line arguments
int process_args(int argc, char* argv[], int* t_eq, int* seed_init, int* seed_bolt, int* seed_spin, double* Tmin, double* Tmax, double* dT) {
  if (argc ==1) {
    //If no parameters passed from command line then set defaults
    *t_eq = 0;
    *seed_init = 100;
    *seed_bolt = 100;
    *seed_spin = 100;
    *Tmin = 1.0;
    *Tmax = 3.4;
    *dT = 0.1;
    std::cout<<"#Using default arguments: "<<*t_eq<<", "<<*seed_init<<", "<<*seed_bolt<<", "<<*seed_spin<<", "<<*Tmin<<", "<<*Tmax<<", "<<dT<<"\n";
  }

  if (argc >1 && argc !=8) {
    //If the number of arguments is awkward, give up and exit the program
    std::cerr<< "\n Invalid set of arguments. Usage: ./ising t_eq seed_init seed_bolt seed_spin T\n";
    exit(1);
  }

  if (argc !=1) {
    //With the correct number of args, put into variables.
    int i = sscanf(argv[1],"%d",t_eq);
    i += sscanf(argv[2],"%d",seed_init);
    i += sscanf(argv[3],"%d",seed_bolt);
    i += sscanf(argv[4],"%d",seed_spin);
    i += sscanf(argv[5],"%la",Tmin);
    i += sscanf(argv[6],"%la",Tmax);
    i += sscanf(argv[7],"%la",dT);
    if (i!=7) {
    	//Unless all scans were successful, throw a tantrum and exit
      std::cerr<<"Failed to fill variables with command line arguments\n Program will now exit\n";
      exit(1);
    }
    //If successful, output the arguments used
    std::cout<<"#Using arguments: "<<*t_eq<<", "<<*seed_init<<", "<<*seed_bolt<<", "<<*seed_spin<<", "<<*Tmin<<", "<<*Tmax<<", "<<dT<<"\n";
  }
  return 0;
}



void output_spin_img(int spins[Ni][Nj]){
  for(int i=0; i<Ni; i++){
    std::cout<<" .";
    for(int j=0; j<Nj; j++){
      if(spins[i][j]==1){std::cout<<"O ";}
      if(spins[i][j]==-1){std::cout<<"+ ";}
      if(spins[i][j]!=-1 && spins[i][j]!=1){std::cout<<"6 ";}
    }
    std::cout<<". \n";
  }
  std::cout<<"\n";
}



void output_parameters(int seed_init, int seed_bolt, int seed_spin, int t_eq){
  std::cout<<"#System params\n";
  std::cout<<"#Seeds (init, bolt, spin):\t("<<seed_init<<", "<<seed_bolt<<", "<<seed_spin<<")\n";
  std::cout<<"#Dimensions (x,y):\t("<<Ni<<", "<<Nj<<")\n";
  std::cout<<"#Properties:\tJ="<<J<<", mu="<<mu<<", H="<<H<<"\n";
  std::cout<<"#Tmax="<<tmax<<",\tTeq="<<t_eq<<"\n";
}


#endif
