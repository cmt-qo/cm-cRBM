//Author: Agnes Valenti. 
//Parts of code provided by G. Carleo and M. Troyer have been utilized, 
//see following copyright notice:

/*
############################ COPYRIGHT NOTICE ##################################

Code provided by G. Carleo and M. Troyer, written by G. Carleo, December 2016.

Permission is granted for anyone to copy, use, modify, or distribute the
accompanying programs and documents for any purpose, provided this copyright
notice is retained and prominently displayed, along with a complete citation of
the published version of the paper:
 ______________________________________________________________________________
| G. Carleo, and M. Troyer                                                     |
| Solving the quantum many-body problem with artificial neural-networks        |
|______________________________________________________________________________|

The programs and documents are distributed without any warranty, express or
implied.

These programs were written for research purposes only, and are meant to
demonstrate and reproduce the main results obtained in the paper.

All use of these programs is entirely at the user's own risk.

################################################################################
*/


#include <vector>
#include <random>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sys/time.h>
//#include <eigen3/Eigen/Dense>
//#include <eigen3/Eigen/Sparse>
//#include <eigen3/Eigen/Core>

#include "cRBM.cc"
#include "AsBp.cc"

//Simple Monte Carlo sampling of a cRBM
const int nsweeps=50000;

//cRBM parameters: a b W ac Wb ap Wp bl Wl Wlc ad
const int npar=2+2+18*2+4+18*2*2+1+9*2+2+3*2+3*4+2;   

template<class Wf,class Hamiltonian> class Sampler{

  //wave-function
  Wf & wf_;

  //Hamiltonian
  Hamiltonian & hamiltonian_;

  //number of spins
  const int nspins_;
  int nv_;
  int nh_;
  int alpha_;
  int npar_;
  int L_;
  //current state in the sampling
  std::vector<int> state_;

  //random number generators and distributions
  std::mt19937 gen_;
  std::uniform_real_distribution<> distu_;
  std::uniform_int_distribution<> distn_;
  std::uniform_int_distribution<> distl_;
  std::uniform_int_distribution<> distp_;
  std::uniform_int_distribution<> distp2_;


  //container for indices of randomly chosen spins to be flipped
  std::vector<int> flips_;

  //useful variable to specify MC updates 
  int ncloop_;

  //quantities needed by the hamiltonian
  //non-zero matrix elements
  std::vector<std::complex<double> > mel_;

  //flip connectors for the hamiltonian(see below for details)
  std::vector<std::vector<int> > flipsh_;

  
  //storage for measured values of the energy
  std::complex<double> en2;
  std::vector<std::complex<double> > Eloc;

  
public:

  Sampler(Wf & wf,Hamiltonian & hamiltonian,int seed):
          wf_(wf),hamiltonian_(hamiltonian),distu_(0,1),nspins_(wf.Nspins()), L_(hamiltonian.getL()),distn_(0,nspins_-1),distl_(0,nspins_/2-1),distp_(0,1),distp2_(0,L_-1)
  {

    Seed(seed);
    nv_=wf_.Nspins();
    nh_=wf_.Nhidden();
    alpha_=wf.alpha();
    npar_=npar;

    Eloc.resize(nsweeps);

  }

 
  //Uniform random number in [0,1)
  inline double Uniform(){
    return distu_(gen_);
  }

  //seed random number generator
  inline void Seed(int seed){
    if(seed<0){
      gen_.seed(std::time(nullptr));
    }
    else{
      gen_.seed(seed);
    }
  }

  //return non-contractible x-loop
  void get_xloop(int d, int number){
    if (d==0){
      for (int i=0; i<L_; i++){
        flips_[i]=number*L_*2+2*i;}
    }
    else if (d==1){
      for (int i=0; i<L_; i++){
        flips_[i]=2*number+1+2*L_*i;}
    }
  }

  //Random spin flips
  //contains single spin updates, vertex flips and non-contractible x-loop flips
  inline bool RandSpin(){
    
    int randnumber=distp_(gen_);
    
    //single spin flip
    if (randnumber==0){
      flips_.resize(1);
      flips_[0]=distn_(gen_);
      ncloop_=0;
    }
    
    else{
      int r2=distp_(gen_);

      //non-contractible loop flip
      if (r2==0){
        flips_.resize(L_);
        int r3=distp2_(gen_);
        get_xloop(distp_(gen_),r3);
        ncloop_=1;
      }

      //vertex flip
      else if(r2==1){
        
        flips_.resize(4);
        int l,o,r,u;  
        int vertex=distl_(gen_);
        hamiltonian_.get_indices_star(vertex,l,o,r,u);
        flips_[0]=l;
        flips_[1]=o;
        flips_[2]=r;
        flips_[3]=u;
        ncloop_=0;
      }
    }    
    return true;
  }

  //Initializes a random spin state
  void InitRandomState(bool mag0=false){
    state_.resize(nspins_);
    for(int i=0;i<nspins_;i++){
      state_[i]=(Uniform()<0.5)?(-1):(1);
    }
  }

 
  //Initializes a 'cold' spin state with all spins up
  void InitColdState(){
    state_.resize(nspins_);
    for(int i=0;i<nspins_;i++){
      state_[i]=1;
    }
  }


  //perform one Monte Carlo step
  void Move(){

    //Picking "nflips" random spins to be flipped
    if(RandSpin()){

      //Computing acceptance probability
      double acceptance=std::norm(wf_.PoP(state_,flips_,ncloop_));

      //Metropolis-Hastings test
      if(acceptance>Uniform()){

        //Updating look-up tables in the wave-function
        wf_.UpdateLt(state_,flips_,ncloop_);

        //Moving to the new configuration
        for(const auto& flip : flips_){
          state_[flip]*=-1;
        }
      }
    }
  }


  //run Monte Carlo sampling to compute energies
  void Run_MC(double thermfactor=0.07,int sweepfactor=1){

    //initialize spin-configuration in cold state
    InitColdState();
    
    //initializing look-up tables in the wave-function
    wf_.InitLt(state_); 
    
    //thermalization
    for(double n=0;n<nsweeps*thermfactor;n+=1){
      for(int i=0;i<nspins_*sweepfactor;i++){  
        Move();
      }
    }
 
  
    //sequence of sweeps
    for(double n=0;n<nsweeps;n+=1){
      for(int i=0;i<nspins_*sweepfactor;i++){
        Move();
      }
      
      //measure local energy of current spin-configuration
      //local energy is saved in member variable en2
      MeasureEnergy();
     
      Eloc[n]=en2;
     
    }
  }


  //save local energy of current spin-configuration in member variable en2
  void MeasureEnergy(){
    en2=0.;

    //Finds the non-zero matrix elements of the hamiltonian
    //on the given state
    //i.e. all the state' such that <state'|H|state> = mel(state') \neq 0
    //state' is encoded as the sequence of spin flips to be performed on state
    hamiltonian_.FindConn(state_,flipsh_,mel_);
    
    for(int i=0;i<flipsh_.size();i++){
      en2+=wf_.PoP(state_,flipsh_[i],0)*mel_[i];
    }
  }
  

 
  //calculate energy and MC sampling error
  void OutputEnergy(double& energy, double& var){
    int nblocks=100;

    int blocksize=std::floor(double(Eloc.size())/double(nblocks));

    double enmean=0;
    double enmeansq=0;

    double enmean_unblocked=0;
    double enmeansq_unblocked=0;

    for(int i=0;i<nblocks;i++){
      double eblock=0;
      for(int j=i*blocksize;j<(i+1)*blocksize;j++){
        eblock+=Eloc[j].real();
        assert(j<Eloc.size());

        double delta=Eloc[j].real()-enmean_unblocked;
        enmean_unblocked+=delta/double(j+1);
        double delta2=Eloc[j].real()-enmean_unblocked;
        enmeansq_unblocked+=delta*delta2;
      }
      eblock/=double(blocksize);
      double delta=eblock-enmean;
      enmean+=delta/double(i+1);
      double delta2=eblock-enmean;
      enmeansq+=delta*delta2;
    }

    enmeansq/=(double(nblocks-1));
    enmeansq_unblocked/=(double((nblocks*blocksize-1)));

    double estav=enmean;
    double esterror=std::sqrt(enmeansq/double(nblocks));

    int ndigits=std::log10(esterror);
    if(ndigits<0){
      ndigits=-ndigits+2;
    }
    else{
      ndigits=0;
    }
    energy=estav;
    var=esterror;
    double t_autocorr=0.5*double(blocksize)*enmeansq/enmeansq_unblocked;
 
  } 

};




