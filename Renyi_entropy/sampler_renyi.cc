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


//Calculate second Renyi entropy via Monte Carlo sampling 
//of a spin Wave-Function

//number of MC sweeps
const int nsweeps=80000;


template<class Wf> class Sampler_Renyi{

  //wave-function (2 copies)
  Wf & wf1_;
  Wf & wf2_;


  //number of spins
  const int nspins_;

  //lattice size: nspins_=2*L_*L_
  int L_;

  //determines size of bipartition
  int A_;

  //useful parameter to track, which copy of the wave function is updated
  int statenum_;

  //current state in the sampling: two copies of the lattice + versions with swapped bipartions
  std::vector<int> state1_;
  std::vector<int> state2_;
  std::vector<int> statemix1_;
  std::vector<int> statemix2_;


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
  
  //storage for measured values of the entropy
  std::complex<double> en2;
  std::vector<std::complex<double> > Eloc;

  
public:

  Sampler_Renyi(Wf & wf1, Wf & wf2, int L,int seed):
          wf1_(wf1),wf2_(wf2),distu_(0,1),nspins_(wf1.Nspins()), L_(L),distn_(0,nspins_-1),distl_(0,nspins_/2-1),distp_(0,1),distp2_(0,L_-1)
  {

   
    Seed(seed);
    
    //set size of the bipartition
    A_=L_*L_;

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

  //returns index of kronecker product
  inline int get_index1D(int x,int y,int AB){
    int z;
    z=L_*2*x+y*2+AB;
    return z;
  }

  //returns spin indices around a vertex v into the variables l1D, o1D, r1D, u1D
  inline void get_indices_star(int v,int &l1D,int &o1D, int &r1D,int &u1D){
    int x,y;
    x=int(v/(L_));
    y=v-L_*x;
    int ym1,xp1;
    ym1=y-1;
    xp1=x+1;
    if (y==0){
       ym1=L_-1;}
    if (x==(L_-1)){
       xp1=0;} 
    l1D=get_index1D(x,ym1,1);
    o1D=get_index1D(x,y,0);
    r1D=get_index1D(x,y,1);
    u1D=get_index1D(xp1,y,0);
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

  //Random spin flips on two copies of state
  //contains single spin updates, vertex flips and non-contractible x-loop flips
  inline bool RandSpin(){
    int randnumber=distp_(gen_);
    int randnumber2=distp_(gen_);

    //update first copy of state
    if (randnumber2==0){
      statenum_=1;
    }

    //update second copy of state
    else{
      statenum_=2;
    }

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
        //std::cout<<"r2 xloop ";
        ncloop_=1;
      }

      //vertex flip
      else if(r2==1){
        
        flips_.resize(4);
        int l,o,r,u;  
        int vertex=distl_(gen_);
        get_indices_star(vertex,l,o,r,u);
        flips_[0]=l;
        flips_[1]=o;
        flips_[2]=r;
        flips_[3]=u;
        ncloop_=0;
      }
    }    
    return true;
  }

  
  //Initializes the first copy as 'cold' spin state with all spins up
  void InitColdState1(){
    state1_.resize(nspins_);
    for(int i=0;i<nspins_;i++){
      state1_[i]=1;
    }
  }

  //Initializes the second copy as 'cold' spin state with all spins up
  void InitColdState2(){
    state2_.resize(nspins_);
    for(int i=0;i<nspins_;i++){
      state2_[i]=1;
    }
  }

  //Initializes the first copy with swapped bipartitions as 'cold' spin state with all spins up
  void InitColdState1mix(){
    statemix1_.resize(nspins_);
    for(int i=0;i<nspins_;i++){
      statemix1_[i]=1;
    }
  }

  //Initializes the second copy with swapped bipartitions as 'cold' spin state with all spins up
  void InitColdState2mix(){
    statemix2_.resize(nspins_);
    for(int i=0;i<nspins_;i++){
      statemix2_[i]=1;
    }
  }


  //perform one Monte Carlo step
  void Move(){

    //Picking "nflips" random spins to be flipped
    if(RandSpin()){

      //Computing acceptance probability
      double acceptance=0;
      if (statenum_==1){
        acceptance=std::norm(wf1_.PoP(state1_,flips_,ncloop_));
      }
      else{
        acceptance=std::norm(wf2_.PoP(state2_,flips_,ncloop_));
      }
          
      //Metropolis-Hastings test
      if(acceptance>Uniform()){

        //Updating look-up tables in the wave-function
        if (statenum_==1){
          wf1_.UpdateLt(state1_,flips_,ncloop_);
        }
        else{
          wf2_.UpdateLt(state2_,flips_,ncloop_);
        }

        //Moving to the new configuration
        if (statenum_==1){
          for(const auto& flip : flips_){
            state1_[flip]*=-1;
            if (flip<A_){
              statemix1_[flip]*=-1;
            }
            else{
              statemix2_[flip]*=-1;
            }
          }
        }
        else{
          for(const auto& flip : flips_){
            state2_[flip]*=-1;
            if (flip<A_){
              statemix2_[flip]*=-1;
            }
            else{
              statemix1_[flip]*=-1;
            }
          }
        }
      }
    }
  }

  

  //run Monte Carlo sampling to compute energies
  void Run_MC(double thermfactor=0.07,int sweepfactor=1,int nflipss=-1){

    //initialize spin-configuration in cold state
    InitColdState1();
    InitColdState2();
    InitColdState1mix();
    InitColdState2mix();
   
    //initializing look-up tables in the wave-function
    wf1_.InitLt(state1_);
    wf2_.InitLt(state2_);
    
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
         
      //measure 'local' entropy of current spin-configuration
      //local entropy is saved in member variable en2  
      MeasureEntropy();
     
      Eloc[n]=en2;
     
      }
    }


  //save local entropy of current spin-configuration in member variable en2
  void MeasureEntropy(){
    en2=0.;

    en2=std::exp(wf1_.LogVal(statemix1_))*std::exp(wf1_.LogVal(statemix2_))/(std::exp(wf1_.LogVal(state1_))*std::exp(wf2_.LogVal(state2_)));
  }
 


  //calculate entropy and MC sampling error
  void OutputEntropy(double& energy, double& var){
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



