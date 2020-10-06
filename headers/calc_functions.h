#ifndef CALC_FUNCTIONS_H
#define CALC_FUNCTIONS_H

typedef double QUADFUNC(double, void*); //for calculating analytic E(T) by gsl integration

//Function declaration and description
//====================================
//functions to calculate quantities
double calc_total_E(const int spins[Ni][Nj]); //uses spin state and global variables to calculate the total energy of the system
double calc_dE(const int spins[Ni][Nj], const int i, const int j); //calculates change in system energy if i,j spin is flipped
int calc_total_S(const int spins[Ni][Nj]); //calculates av magnetisation
void calc_av_S(int S[tmax], double &avS, double &devS); //calculate average S and the std error in it
void calc_av_E(double E[tmax], double &avE, double &devE, int t_eq);
void calc_dEdT(double E1, double E2, double dT, double errE, double &C, double &errC);
double calc_exact_M_perspin(double T); //returns Onsager's analytical value (to computer precision) for given T
double calc_exact_E_perspin(double T, double* err); 
double calc_exact_c_perspin(double T, double* err);
double est_exact_c_perspin(double T);

//bits and pieces used to calculate the analytic values of E and C
double K(double T); //parameter used in exact E (strictly at H=0 field)
double dKd1oT(double T); //first deriv of K wrt 1/T
double d2Kd1oT(double T); //second deriv
double integrand(double theta, void* T);
double diffd_integrand(double x, void* params);
double integrate( QUADFUNC func, void* params, double x0, double x1, double & error);

//Function implementations
//====================================
inline double calc_total_E(const int spins[Ni][Nj]){
  int interaction = 0;
  int ext = 0;
  int imax=0;
  int jmax=0;
  for(int i=0; i<Ni; i++){
    for(int j=0; j<Nj; j++){
      if(i!=Ni-1){imax=i+1;}
      else{imax=0;}
      if(j!=Nj-1){jmax=j+1;}
      else{jmax=0;}
      
      interaction += spins[i][j]*(spins[imax][j] + spins[i][jmax]);
      ext += spins[i][j];
    }
  }
    
  return ( -((double)J)*((double)interaction) - mu*H*((double)ext) );
}



inline double calc_dE(const int spins[Ni][Nj], const int i, const int j){  
  //Copy the spin under scrutiny
  //int test_spin = spins[i][j];
  
  //Find interaction energy before  
  unsigned int imin;
  unsigned int jmin;
  unsigned int imax;
  unsigned int jmax;
  if(i!=0){imin=i-1;} //% operator leaves negative values negative
  else{imin=Ni-1;}
  if(j!=0){jmin=j-1;}
  else{jmin=Nj-1;}

  if(i!=Ni-1){imax=i+1;}
  else{imax=0;}
  if(j!=Nj-1){jmax=j+1;}
  else{jmax=0;}

  int startEi = spins[i][j]*(spins[imin][j] + spins[imax][j] + spins[i][jmin] + spins[i][jmax]);

  //Now calculate energy associated with external field interaction (remember spin already flipped)
  //int dSi = -2*test_spin;
  
  //double dEh = -mu*H*((double)dSi);
  
  //So the total change is
  return (J*((double)(2*startEi)));//(dEi+dEh);

}



inline int calc_total_S(const int spins[Ni][Nj]){
  int S=0;
  for (int i=0;i<Ni;i++){
    for (int j=0;j<Nj;j++){
      S =S+ spins[i][j];
    }
  }
  return S;
}



void calc_av_S(int S[tmax], double &avS, double &devS, int t_eq){
  avS=0;
  devS=0;
  double sum=0;
  double sumsq=0;
  for(int i=t_eq; i<tmax; i++){
    sum+=S[i];
    sumsq+=S[i]*S[i];
  }
  avS=((double)sum/((double)(tmax-t_eq)));
  double avSsq = ((double)sumsq/((double)(tmax-t_eq)));
  devS = sqrt(avSsq-avS*avS);
}



void calc_av_E(double E[tmax], double &avE, double &devE, int t_eq){
  avE=0;
  devE=0;
  double sum=0;
  double sumsq=0;
  for(int i=t_eq; i<tmax; i++){
    sum+=E[i];
    sumsq+=E[i]*E[i];
  }
  avE=((double)sum/((double)(tmax-t_eq)));
  double avEsq = ((double)sumsq/((double)(tmax-t_eq)));
  devE = sqrt(avEsq-avE*avE);
}



double calc_exact_M_perspin(double T){
  //M(T) = (1-[sinh(2J/kbT)]^-4)^1/8
  double x=2.0/T;
  return pow((1 - pow(sinh(x),-4.0) ), 0.125);
}



double K(double T){
  //K=2/(cosh(2J/kbT)coth(2J/kbT))
  return ( (2.0*tanh(2.0/T))/cosh(2.0/T) );
}



double dKd1oT(double T){
  //dKd1oT = 4*(1-sinh(2/T)^2)/cosh(2/T)^3
  return ( 4*(1-sinh(2.0/T)*sinh(2.0/T))/(pow(cosh(2.0/T), 3)) );
}



double d2Kd1oT(double T){
  //second deriv of K wrt 1/T = 4*(cosh(4/T)-11)sinh(2/T)/cosh^4(2/T)
  return ( 4.0*(cosh(4.0/T)-11)*sinh(2.0/T)/pow(cosh(2.0/T),4.0) );
}

double integrand(double x, void* params){
  double T = *(double*)params;
  double numerator = (sin(x)*sin(x));
  double D = sqrt( 1.0-pow((K(T)*sin(x)),2.0) );
  return (numerator/(D*(1.0+D)));
}

double diffd_integrand(double x, void* params){
  //integrand: sin^2(x)/( sqrt(1-K^2sin^2(x))*(1+sqrt(1-K^2sin^2(x)) = sin^2(x)/D(1+D)
  //differentiated wrt to 1/T : -sin^2(x)(1+2D/(D(1+D))^2)*dD/d(1/T) = A*dDd1oT
  //dD/d(1/T) = -(Ksin^2(x)/D)*dKd1oT
  double T = *(double*)params;
  double D = sqrt( 1.0-pow((K(T)*sin(x)),2.0) );
  double A = ( -sin(x)*sin(x)*(1.0+2.0*D) )/(pow((D*(1.0+D)),2));
  double dDd1oT = -( (K(T)*sin(x)*sin(x))/D )*dKd1oT(T);
  return (A*dDd1oT);
	
}

double integrate( QUADFUNC func, void* params, double x0, double x1, double & error){
  const int max_intervals = 1000000;
  gsl_integration_workspace* workspace = gsl_integration_workspace_alloc( max_intervals);
  double result;
  gsl_function F;
  F.function = func;
  F.params = params;
  gsl_integration_qag( &F, x0, x1, 0, 1e-7, max_intervals, 6, workspace, &result, &error);
  gsl_integration_workspace_free(workspace);
  return result;
}

double calc_exact_E_perspin(double T, double* err){
  //E=-2Jtanh(2J/kbT) + K/2pi * dK/d(1/kb*T) * integral[0,pi] of dx*sin^2(x)/( sqrt(1-K^2sin^2(x))*(1+sqrt(1-K^2sin^2(x)) )
  //Will do some numerical integration
  double meanfield = -2.0*J*tanh(2.0/T);
  double pi = C::pi;
  double prefactor = (K(T)/(2.0*pi)) * dKd1oT(T);
	
  double error=0;
  double integral = integrate(integrand, &T, 0, pi, error);
	
  *err = error;
  return (meanfield + prefactor*integral);
}

double calc_exact_c_perspin(double T, double* err){
  //C=dE/dT
  //C = (4J/T^2)*sech^2(2/T)  +  (-1/2*pi*T^2)*[ (dK/d(1/T))^2 * integral +....
  //+... K* d2K/d(1/T)^2 * integral + KdK/d(1/T) * diffd_integral]
  double diffd_meanfield =  (4.0*J)/(T*T*cosh(2.0/T)*cosh(2.0/T));
  double pi = C::pi;
  double error = 0;
  double integral = integrate(integrand, &T, 0, pi, error);
  double diffd_integral = integrate(diffd_integrand, &T, 0, pi, error);
  *err = error;
	
  double prefactor = (-1.0/(2*pi*T*T));
  double term1 = pow(dKd1oT(T),2.0)*integral;
  double term2 = K(T)*d2Kd1oT(T)*integral;
  double term3 = K(T)*dKd1oT(T)*diffd_integral;
  return ( (diffd_meanfield) + prefactor*(term1 + term2 + term3));	
}

#endif
