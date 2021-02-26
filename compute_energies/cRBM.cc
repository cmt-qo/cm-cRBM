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


#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <complex>
#include <fstream>
#include <cassert>
#include "useful_functions.cc"
#include<algorithm> // for copy() and assign() 
#include<iterator> // for back_inserter 

//cRBM constructed for the disordered toric model
class cRBM{

  //Neural-network weights connected to single spin visible neurons
  std::vector<std::vector<std::complex<double> > > W_;

  //Neural-network weights connected to bond visible neurons
  std::vector<std::vector<std::complex<double> > > Wb_;

  //Neural-network weights connected to plaquette visible neurons
  std::vector<std::vector<std::complex<double> > > Wp_;

  //Neural-network weights connected to non-contractible loop correlators and extra hidden neurons
  std::vector<std::vector<std::complex<double> > > Wl_;

  //Neural-network weights connected to non-contractible loop correlators
  std::vector<std::vector<std::complex<double> > > Wlc_;

  //Neural-network visible bias (single spin visible neurons)
  std::vector<std::complex<double> > a_;

  //Neural-network visible bias (bond correlators)
  std::vector<std::complex<double> > ac_;
  
  //Neural-network visible bias (plaquette correlators)
  std::vector<std::complex<double> > ap_;
  
  //Neural-network visible bias (non-contractible loop correlators)
  std::vector<std::complex<double> > ad_;

  //Neural-network hidden bias
  std::vector<std::vector<std::complex<double> > > b_;

  //Neural-network hidden bias extra hidden neurons, 
  //which are only connected to non-contractible loop correlators
  std::vector<std::vector<std::complex<double> > > bl_;

 
  //Number of hidden units*number of symmetries
  int nh_;

  //Number of hidden units
  int alpha_;

  //Number of spins
  int nv_;

  //lattice length, nv_=2*L*L
  int L;

  //look-up tables
  std::vector<std::complex<double> > Lt_;

  //look-up tables for extra hidden neurons
  std::vector<std::complex<double> > Ltl_;

  //Useful quantities for safe computation of ln(cosh(x))
  const double log2_;

public:

  //vectors to store symmetry maps
  std::vector<std::vector<int> > symm;
  std::vector<std::vector<int> > symm_inv;

  cRBM(std::string filename, std::string filename_symm):log2_(std::log(2.)){
    LoadParameters(filename);
    Init_symm(filename_symm);
  }
  

//save spin indices of the two spins connected to bond b_ind into the variables i1,i2
inline void get_bond_indices(int b_ind, int &i1, int &i2){
  //x coordinate 
  int x=b_ind/(2*2*L);

  //y coordinate
  int y=(b_ind-x*2*2*L)/(2*2);

  //AB sublattice
  int AB=(b_ind-x*2*2*L-4*y)/2;

  //bond sublattice
  int b=b_ind-x*2*2*L-4*y-2*AB;

  i1=2*L*x+2*y+AB;
  i2=2*L*((x-1+AB+L)%L)+2*((y-1+AB+b+L)%L)+(1-AB);
}

//compute adjacent bond indices (and spin indices) b1,b2,b3,b4 (i1,i2,i3,i4) 
//connecting the spin 'flip' to the spins i1,i2,i3,i4
inline void get_bonds(int flip,int &i1,int &i2,int &i3,int &i4,int &b1,int &b2,int &b3,int &b4){
  int x=flip/(2*L);
  int y=(flip-x*2*L)/(2);
  int AB=(flip-x*2*L-2*y);
  b1=4*L*x+4*y+2*AB;
  b2=b1+1;
  b3=4*L*((x+AB+L)%L)+4*((y-1+AB+L)%L)+(1-AB)*2+1;
  b4=4*L*((x+AB+L)%L)+4*((y+AB+L)%L)+(1-AB)*2;
  i1=2*L*((x-1+AB+L)%L)+2*((y-1+AB+L)%L)+(1-AB);
  i2=2*L*((x-1+AB+L)%L)+2*((y+AB+L)%L)+(1-AB);
  i3=2*L*((x+AB+L)%L)+2*((y-1+AB+L)%L)+(1-AB);
  i4=2*L*((x+AB+L)%L)+2*((y+AB+L)%L)+(1-AB);
}

//return z-loop correlator 'j' s_1*...*s_k shifted by symmetry operation 's'
inline int get_zloop(const std::vector<int> & state, int j, int s){
  int loop=1;
  if (s==-1){
    if (j<L){   
      for (int i=0; i<L; i++){
        loop=loop*state[2*L*i+2*j];
      }
    }
    else{    
      for (int i=0; i<L; i++){
        loop=loop*state[2*i+1+2*L*(j-L)];
      }
    }
  }

  else{
    if (j<L){   
      for (int i=0; i<L; i++){
        loop=loop*state[symm[2*L*i+2*j][s]];
      }
    }
    else{    
      for (int i=0; i<L; i++){
        loop=loop*state[symm[2*i+1+2*L*(j-L)][s]];
      }
    }
  }

  return loop;
} 

//return z-loop index, of which the spin 'spin' is part of
inline int spin_to_zloop(int spin){
  int z=-1;
  int x=int(spin/(2*L));
  int y=int((spin-2*L*x)/2);
  int AB=spin-2*L*x-2*y;
  if (AB==0)
    z=y;
  else if (AB==1)
    z=L+x;
  return z;
}


//shift index of non-contractible loop i by symmetry operator s
inline int symm_L(int i, int s){
  return (i+s)%L;
}

//shift index of non-contractible loop i by inverted action of symmetry operator s
inline int symm_inv_L(int i, int s){
  return (i-s+L)%L;
}

//computes the logarithm of the wave-function
inline std::complex<double> LogVal(const std::vector<int> & state){

  std::complex<double> rbm(0.,0.);

  for(int v=0;v<nv_;v++){
    rbm+=a_[v%2]*double(state[v]);
  }

  for(int vb=0;vb<nv_*2;vb++){
    int i1_,i2_;
    get_bond_indices(vb, i1_,i2_);
    rbm+=ac_[vb%4]*double(state[i1_])*double(state[i2_]);
  }

  for(int vp=0;vp<nv_/2;vp++){
    int l,o,r,u;
    get_indices_plaquette_(vp,l,o,r,u,L);
    rbm+=ap_[0]*double(state[l])*double(state[o])*double(state[r])*double(state[u]);
  }

   
  for(int vd=0;vd<2*L;vd++){
    if (vd<2*L){
      rbm+=ad_[int(vd/L)]*double(get_zloop(state, vd, -1));
    }  
  }


  for(int ai=0; ai<alpha_; ai++){  
    for(int s=0;s<nv_/2;s++){    //hier nv_ in Rolle von nh_/alpha
      std::complex<double> thetah=b_[0][ai];
      for(int v=0;v<nv_;v++){
        thetah+=double(state[symm[v][s]])*(W_[v][ai]);
      }
      
      for(int bi=0; bi<nv_*2; bi++){
        int i1_,i2_;
        get_bond_indices(bi, i1_,i2_);
        thetah+=double(state[symm[i1_][s]])*double(state[symm[i2_][s]])*(Wb_[bi][ai]);
      }

      for(int pi=0; pi<nv_/2; pi++){
        int l,o,r,u;
        get_indices_plaquette_(pi,l,o,r,u,L);
        thetah+=double(state[symm[l][s]])*double(state[symm[o][s]])*double(state[symm[r][s]])*double(state[symm[u][s]])*(Wp_[pi][ai]);
      }
      for(int l=0;l<2*L;l++){
        double loopvalue=double(get_zloop(state, l,s));
        thetah+=loopvalue*(Wlc_[l][ai]);
      }
      rbm+=cRBM::lncosh(thetah);
    }
  }

  for (int ai=0; ai<2; ai++){
    for(int s=0; s<L; s++){
      std::complex<double> thetah=bl_[0][ai];
      for(int l=0;l<L;l++){
        double loopvalue=double(get_zloop(state, ai*L+symm_L(l,s),-1));
        thetah+=loopvalue*(Wl_[l][ai]);
        //std::cout<<"ai: "<<ai<<" s: "<<s<<" l: "<<l<<" loop insges: "<<ai*L+symm_L(l,s)<< "loop get zloop: "<<loopvalue<<std::endl;
      }
      rbm+=cRBM::lncosh(thetah);
    }
  }
  return rbm;
}

//computes the logarithm of Psi(state')/Psi(state)
//where state' is a state with a certain number of flipped spins
//the vector "flips" contains the sites to be flipped
//look-up tables are used to speed-up the calculation
//if the spin-flip corresponds to a non-contractible x-loop, ncloop=1
inline std::complex<double> LogPoP(const std::vector<int> & state, const std::vector<int> & flips, const int ncloop){

  if(flips.size()==0){
    return 0.;
  }

  std::complex<double> logpop(0.,0.);

  //Change due to the visible bias on single-spin visible neurons
  for(const auto & flip : flips){
    logpop-=a_[flip%2]*2.*double(state[flip]);
  }

  //Change due to the visible bias on non-contractible z-loop correlators
  if (flips.size()==1 || ncloop==1){
    for(const auto & flip : flips){
      int ci=spin_to_zloop(flip);
      int loop=get_zloop(state,ci,-1);    
      logpop-=2.*ad_[int(ci/L)]*double(loop); 
    }
  }

  //single-spin flips (distinguish updates for different types of spin flips)
  if (flips.size()==1){

    //Change due to the visible bias on bond correlators
    for(const auto & flip : flips){
      int i1,i2,i3,i4,b1,b2,b3,b4;
      get_bonds(flip,i1,i2,i3,i4,b1,b2,b3,b4);
      logpop-=ac_[b1%4]*2.*double(state[flip])*double(state[i1]);
      logpop-=ac_[b2%4]*2.*double(state[flip])*double(state[i2]);
      logpop-=ac_[b3%4]*2.*double(state[flip])*double(state[i3]);
      logpop-=ac_[b4%4]*2.*double(state[flip])*double(state[i4]);
    }

    //Change due to the visible bias on plaquette correlators
    for(const auto & flip : flips){
      int p1,p2,l1,o1,r1,u1,l2,o2,r2,u2;
      get_2plaquettes(flip,p1,p2,L);
      get_indices_plaquette_(p1,l1,o1,r1,u1,L);
      get_indices_plaquette_(p2,l2,o2,r2,u2,L);
      logpop-=ap_[0]* 2.*double(state[l1])*double(state[o1])*double(state[r1])*double(state[u1]);
      logpop-=ap_[0]* 2.*double(state[l2])*double(state[o2])*double(state[r2])*double(state[u2]);
    }
  }

  //vertex flips or non-contractible loop flips (distinguish updates for different types of spin flips)
  else if (flips.size()>1){

    //Change due to the visible bias on bond correlators
    std::vector<int> bonds(4*L*L,1);
    for(const auto & flip : flips){
      int i1,i2,i3,i4,b1,b2,b3,b4;
      get_bonds(flip,i1,i2,i3,i4,b1,b2,b3,b4);
      logpop-=ac_[b1%4]*2.*double(state[flip])*double(state[i1])*double(bonds[b1]);
      bonds[b1]=-1;
      logpop-=ac_[b2%4]*2.*double(state[flip])*double(state[i2])*double(bonds[b2]);
      bonds[b2]=-1;
      logpop-=ac_[b3%4]*2.*double(state[flip])*double(state[i3])*double(bonds[b3]);
      bonds[b3]=-1;
      logpop-=ac_[b4%4]*2.*double(state[flip])*double(state[i4])*double(bonds[b4]);
      bonds[b4]=-1;
    }
  }


  //Change due to the interaction weights
  for(int h=0;h<nh_;h++){

    //hidden neuron index
    int ai=int(h/(nv_/2));

    //symmetry index
    int s=h-(nv_/2)*ai;

    //use look-up table
    std::complex<double> thetah=Lt_[h];
    std::complex<double> thetahp=thetah;

    //Change due to the interaction weights connected to single-spin visible neurons
    for(const auto & flip : flips){
      thetahp-=2.*double(state[flip])*(W_[symm_inv[flip][s]][ai]);
    }

    //Change due to the interaction weights connected to non-contractible z-loop correlators
    if (flips.size()==1 || ncloop==1){
      for(const auto & flip : flips){
          int ci=spin_to_zloop(flip);
          int bWloop=spin_to_zloop(symm_inv[flip][s]);
          double loopvalue=double(get_zloop(state, ci,-1));
          thetahp-=2.*loopvalue*(Wlc_[bWloop][ai]);
      }
    }

    //single spin-flips (distinguish updates for different types of spin flips)
    if (flips.size()==1){

      //Change due to the interaction weights connected to bond correlators
      for(const auto & flip : flips){
        int i1,i2,i3,i4,b1,b2,b3,b4,j1,j2,j3,j4;
        get_bonds(flip,i1,i2,i3,i4,b1,b2,b3,b4);
        get_bonds(symm_inv[flip][s],j1,j2,j3,j4,b1,b2,b3,b4);
        thetahp-=Wb_[b1][ai]*2.*double(state[flip])*double(state[i1]);
        thetahp-=Wb_[b2][ai]*2.*double(state[flip])*double(state[i2]);
        thetahp-=Wb_[b3][ai]*2.*double(state[flip])*double(state[i3]);
        thetahp-=Wb_[b4][ai]*2.*double(state[flip])*double(state[i4]);
      }

      //Change due to the interaction weights connected to plaquette correlators
      for(const auto & flip : flips){
        int p1,p2,l1,o1,r1,u1,l2,o2,r2,u2;   
        get_2plaquettes(flip,p1,p2,L);
        get_indices_plaquette_(p1,l1,o1,r1,u1,L);
        get_indices_plaquette_(p2,l2,o2,r2,u2,L);
        get_2plaquettes(symm_inv[flip][s],p1,p2,L);
        thetahp-=Wp_[p1][ai]* 2.*double(state[l1])*double(state[o1])*double(state[r1])*double(state[u1]);
        thetahp-=Wp_[p2][ai]* 2.*double(state[l2])*double(state[o2])*double(state[r2])*double(state[u2]);
      }
    }

    //vertex flips or non-contractible loop flips (distinguish updates for different types of spin flips)
    else if (flips.size()>1){

      //Change due to the interaction weights connected to bond correlators
      std::vector<int> bonds(4*L*L,1);
      for(const auto & flip : flips){
        int i1,i2,i3,i4,b1,b2,b3,b4,j1,j2,j3,j4;
        get_bonds(flip,i1,i2,i3,i4,b1,b2,b3,b4);
        get_bonds(symm_inv[flip][s],j1,j2,j3,j4,b1,b2,b3,b4);
        thetahp-=Wb_[b1][ai]*2.*double(state[flip])*double(state[i1])*double(bonds[b1]);
        bonds[b1]=-1;
        thetahp-=Wb_[b2][ai]*2.*double(state[flip])*double(state[i2])*double(bonds[b2]);
        bonds[b2]=-1;
        thetahp-=Wb_[b3][ai]*2.*double(state[flip])*double(state[i3])*double(bonds[b3]);
        bonds[b3]=-1;
        thetahp-=Wb_[b4][ai]*2.*double(state[flip])*double(state[i4])*double(bonds[b4]);
        bonds[b4]=-1;
      }
    }
    logpop+= ( cRBM::lncosh(thetahp)-cRBM::lncosh(thetah) );
  }
    

  //change due to two extra neurons only connected to non-contractible z-loop correlators
  //no change if vertices are flipped
  if (flips.size()==1 || ncloop==1){
    for(int h=0;h<L*2;h++){

      //hidden neuron index
      int ai=int(h/(L));

      //symmetry index
      int s=h-(L)*ai;

      //look-up table
      std::complex<double> thetah=Ltl_[h];
      std::complex<double> thetahp=thetah;
        
      for(const auto & flip : flips){
        int ci=spin_to_zloop(flip);
        int loop=get_zloop(state,ci,-1);
        if (ci<L && ai==0){
          int b_i=symm_inv_L(ci,s);
          thetahp-=Wl_[b_i][ai]*2.*double(loop);}
        else if (ci>=L && ai==1){
          int b_i=symm_inv_L(ci-L,s);
          thetahp-=Wl_[b_i][ai]*2.*double(loop);}
      } 
      logpop+= ( cRBM::lncosh(thetahp)-cRBM::lncosh(thetah) );  
    }
  }
  return logpop;
}

inline std::complex<double> PoP(const std::vector<int> & state,const std::vector<int> & flips, const int ncloop){
  return std::exp(LogPoP(state,flips,ncloop));
}

//initialization of the look-up tables
void InitLt(const std::vector<int> & state){
  Lt_.resize(nh_);

  for(int h=0;h<nh_;h++){

    //hidden neuron index
    int ai=int(h/(nv_/2));

    //symmetry index
    int s=h-ai*(nv_/2);

    //look-up table: Lt_[h]= hidden bias+sum over all correlators(weights*correlator)
    //hidden bias
    Lt_[h]=b_[0][ai];

    //Contribution due to the interaction weights connected to single-spin visible neurons
    for(int v=0;v<nv_;v++){
      Lt_[h]+=double(state[symm[v][s]])*(W_[v][ai]);
    }

    //Contribution due to the interaction weights connected to bond correlators
    for(int bi=0; bi<nv_*2; bi++){
      int i1_,i2_;
      get_bond_indices(bi, i1_,i2_);
      Lt_[h]+=double(state[symm[i1_][s]])*double(state[symm[i2_][s]])*Wb_[bi][ai];
    }

    //Contribution due to the interaction weights connected to plaquette correlators
    for(int pi=0; pi<nv_/2; pi++){
      int l,o,r,u;
      get_indices_plaquette_(pi,l,o,r,u,L);
      Lt_[h]+=double(state[symm[l][s]])*double(state[symm[o][s]])*double(state[symm[r][s]])*double(state[symm[u][s]])*(Wp_[pi][ai]);
    }

    //Contribution due to the interaction weights connected to non-contractible z-loop correlators
    for(int l=0;l<2*L;l++){
      double loopvalue=double(get_zloop(state, l,s));
      Lt_[h]+=loopvalue*(Wlc_[l][ai]);
    }
      
    }
    

    //look-up tables for two extra hidden neurons
    Ltl_.resize(2*L);
    for (int ai=0; ai<2; ai++){
      for(int s=0; s<L; s++){
        Ltl_[ai*L+s]=bl_[0][ai];
        for(int l=0;l<L;l++){
          double loopvalue=double(get_zloop(state, ai*L+symm_L(l,s),-1));
          Ltl_[ai*L+s]+=loopvalue*(Wl_[l][ai]);
        }
      }
    }
  }

//updates the look-up tables after spin flips
//the vector "flips" contains the indices of sites to be flipped
//if the spin-flip corresponds to a non-contractible x-loop, ncloop=1
void UpdateLt(const std::vector<int> & state,const std::vector<int> & flips, const int ncloop){
    
  if(flips.size()==0){
    return;
  }

    
  for(int h=0;h<nh_;h++){

    //hidden neuron index
    int ai=int(h/(nv_/2));

    //symmetry index
    int s=h-ai*(nv_/2);
      
    //Change due to the interaction weights connected to single-spin visible neurons
    for(const auto & flip : flips){
      Lt_[h]-=2.*double(state[flip])*W_[symm_inv[flip][s]][ai];
    }

    //Change due to the interaction weights connected to non-contractible z-loop correlators
    if (flips.size()==1 || ncloop==1){
      for(const auto & flip : flips){
        int ci=spin_to_zloop(flip);
        int bWloop=spin_to_zloop(symm_inv[flip][s]);
        double loopvalue=double(get_zloop(state, ci,-1));
        Lt_[h]-=2.*loopvalue*(Wlc_[bWloop][ai]);
      }
    }

    //single spin-flips (distinguish updates for different types of spin flips)
    if (flips.size()==1){

      //Change due to the interaction weights connected to bond correlators
      for(const auto & flip : flips){
        int i1,i2,i3,i4,b1,b2,b3,b4,j1,j2,j3,j4;
        get_bonds(flip,i1,i2,i3,i4,b1,b2,b3,b4);
        get_bonds(symm_inv[flip][s],j1,j2,j3,j4,b1,b2,b3,b4);
        Lt_[h]-=Wb_[b1][ai]*2.*double(state[flip])*double(state[i1]);
        Lt_[h]-=Wb_[b2][ai]*2.*double(state[flip])*double(state[i2]);
        Lt_[h]-=Wb_[b3][ai]*2.*double(state[flip])*double(state[i3]);
        Lt_[h]-=Wb_[b4][ai]*2.*double(state[flip])*double(state[i4]);
        }
  
      //Change due to the interaction weights connected to plaquette correlators
      for(const auto & flip : flips){
        int p1,p2,l1,o1,r1,u1,l2,o2,r2,u2;   
        get_2plaquettes(flip,p1,p2,L);
        get_indices_plaquette_(p1,l1,o1,r1,u1,L);
        get_indices_plaquette_(p2,l2,o2,r2,u2,L);
        get_2plaquettes(symm_inv[flip][s],p1,p2,L);
        Lt_[h]-=Wp_[p1][ai]* 2.*double(state[l1])*double(state[o1])*double(state[r1])*double(state[u1]);
        Lt_[h]-=Wp_[p2][ai]* 2.*double(state[l2])*double(state[o2])*double(state[r2])*double(state[u2]);
      }
    }

    //vertex flips or non-contractible loop flips (distinguish updates for different types of spin flips)
    else if (flips.size()>1){
      
      //Change due to the interaction weights connected to bond correlators
      std::vector<int> bonds(4*L*L,1);
      for(const auto & flip : flips){
        int i1,i2,i3,i4,b1,b2,b3,b4,j1,j2,j3,j4;
        get_bonds(flip,i1,i2,i3,i4,b1,b2,b3,b4);
        get_bonds(symm_inv[flip][s],j1,j2,j3,j4,b1,b2,b3,b4);
        Lt_[h]-=Wb_[b1][ai]*2.*double(state[flip])*double(state[i1])*double(bonds[b1]);
        bonds[b1]=-1;
        Lt_[h]-=Wb_[b2][ai]*2.*double(state[flip])*double(state[i2])*double(bonds[b2]);
        bonds[b2]=-1;
        Lt_[h]-=Wb_[b3][ai]*2.*double(state[flip])*double(state[i3])*double(bonds[b3]);
        bonds[b3]=-1;
        Lt_[h]-=Wb_[b4][ai]*2.*double(state[flip])*double(state[i4])*double(bonds[b4]);
        bonds[b4]=-1;
      }
    }
  }


  //change due to two extra neurons only connected to non-contractible z-loop correlators
  //no change if vertices are flipped
  if (flips.size()==1 || ncloop==1){
    for(int h=0;h<L*2;h++){

      int ai=int(h/(L));

      int s=h-(L)*ai;
  
      for(const auto & flip : flips){
        int ci=spin_to_zloop(flip);
        int loop=get_zloop(state,ci,-1);
          
        if (ci<L && ai==0){
          int b_i=symm_inv_L(ci,s);
          Ltl_[h]-=Wl_[b_i][ai]*2.*double(loop);}
          
        else if (ci>=L && ai==1){
          int b_i=symm_inv_L(ci-L,s);
          Ltl_[h]-=Wl_[b_i][ai]*2.*double(loop);}
      } 
    }
  }
}

//read translational symmetry maps from file 'filename_symm' 
//into vectors symm (for action of symmetry operators) and symm_inv (for inverse action of the symmetry operators)
void Init_symm(std::string filename_symm){
  symm.resize(nv_,std::vector<int>(nv_/2));
  symm_inv.resize(nv_,std::vector<int>(nv_/2));
  std::ifstream fin(filename_symm.c_str());

  if(!fin.good()){
    std::cerr<<"# Error : Cannot load from file "<<filename_symm<<" : file not found."<<std::endl;
    std::abort();
  }
  for(int s=0; s<nv_/2; s++){
    for (int i=0; i<nv_; i++){
      fin>>symm[i][s];}
  }

  for(int s=0; s<nv_/2; s++){
    for (int i=0; i<nv_; i++){
      symm_inv[symm[i][s]][s]=i;}
  }
}

//loads the parameters of the cRBM wave-function from a given file
void LoadParameters(std::string filename){

  std::ifstream fin(filename.c_str());

  if(!fin.good()){
    std::cerr<<"# Error : Cannot load from file "<<filename<<" : file not found."<<std::endl;
    std::abort();
  }

  std::complex<double> a_placeholder;
  std::complex<double> b_placeholder;
  std::complex<double> alpha_placeholder;

  fin>>a_placeholder;
  nv_=(int)(a_placeholder.real());

  fin>>b_placeholder;
  nh_=(int)(b_placeholder.real());

  fin>>alpha_placeholder;
  alpha_=(int)(alpha_placeholder.real());

  if(!fin.good() || nv_<0 || nh_<0){
    std::cerr<<"# Trying to load from an invalid file.";
    std::cerr<<std::endl;
    std::abort();
  }

  //lattice length
  L=int(std::sqrt(nv_/2));

  //resize variable vectors
  //visible biases
  a_.resize(2);  
  ac_.resize(4);  
  ap_.resize(1);  
  ad_.resize(2);  

  //hidden biases
  b_.resize(1,std::vector<std::complex<double> > (alpha_)); 
  bl_.resize(1,std::vector<std::complex<double> > (alpha_)); 

  //interaction weights
  W_.resize(nv_,std::vector<std::complex<double> > (alpha_));
  Wb_.resize(nv_*2,std::vector<std::complex<double> > (alpha_));
  Wp_.resize(nv_/2,std::vector<std::complex<double> > (alpha_));
  Wl_.resize(L,std::vector<std::complex<double> > (alpha_));
  Wlc_.resize(2*L,std::vector<std::complex<double> > (alpha_));

  for(int i=0;i<2;i++){
    fin>>a_[i];
  }  

  for(int ai=0; ai<alpha_; ai++){
    for(int j=0;j<1;j++){
      fin>>b_[j][ai];
    }
  }
  
  for(int ai=0;ai<alpha_;ai++){
    for(int i=0;i<nv_;i++){
      fin>>W_[i][ai];
    }
  }
  
  for(int i=0;i<4;i++){
    fin>>ac_[i];
  } 

  for(int ai=0;ai<alpha_;ai++){ 
    for(int i=0;i<nv_*2;i++){
      fin>>Wb_[i][ai];        
    }
  }

  for(int i=0;i<1;i++){
    fin>>ap_[i];
  } 

  for(int ai=0;ai<alpha_;ai++){ 
    for(int i=0;i<nv_/2;i++){
      fin>>Wp_[i][ai];        
    }
  }
  
  for(int ai=0;ai<2;ai++){
    for(int i=0;i<1;i++){
      fin>>bl_[i][ai];  
    } 
  }

  for(int ai=0;ai<2;ai++){ 
    for(int i=0;i<L;i++){
      fin>>Wl_[i][ai];        
    }
  }

  for(int ai=0;ai<2;ai++){ 
    for(int i=0;i<2*L;i++){
      fin>>Wlc_[i][ai];        
    }
  }

  for(int i=0;i<2;i++){
    fin>>ad_[i];
  }  
}

//ln(cos(x)) for real argument
//for large values of x we use the asymptotic expansion
inline double lncosh(double x)const{
  const double xp=std::abs(x);
  if(xp<=12.){
    return std::log(std::cosh(xp));
  }
  else{
    return xp-log2_;
  }
}

inline double lnsinh(double x)const{
  const double xp=std::abs(x);
  if(xp<=12.){
    return std::log(std::sinh(xp));
  }
  else{
    return xp-log2_;
  }
}

//ln(cos(x)) for complex argument
//the modulus is computed by means of the previously defined function
//for real argument
inline std::complex<double> lncosh(std::complex<double> x)const{
  const double xr=x.real();
  const double xi=x.imag();

  std::complex<double> res=cRBM::lncosh(xr);
  res +=std::log( std::complex<double>(std::cos(xi),std::tanh(xr)*std::sin(xi)) );

  return res;
}

//total number of spins
inline int Nspins()const{
  return nv_;
}
 
//return number of hidden neurons*number of symmetry operators
inline int Nhidden()const{
  return nh_;
}

//return number of hidden neurons
inline int alpha()const{
  return alpha_;
}

//return tanh of look-up table Lt_[k]
//k index in hidden layer 0<=k<nh_
inline std::complex<double> get_tanh(int k){   
  return std::tanh(Lt_[k]); 
}

//return tanh of look-up table Ltl_[k]
//k index in hidden layer 0<=k<nh_
inline std::complex<double> get_tanhl(int k){   
  return std::tanh(Ltl_[k]); //?
}

//return argument theta=Lt_[k]
inline std::complex<double> theta(int k){
  return Lt_[k];
}

};



