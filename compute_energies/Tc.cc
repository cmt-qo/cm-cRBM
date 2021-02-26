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
#include <vector>
#include <complex>

//constructs disordered toric Hamiltonian
class Tc{

  //lattice size
  const int L;
 
  //useful variables to save spin indices
  int l,o,r,u;

  //magnetic field strengths
  double h_x,h_y,h_z;

  //stabilizer weight
  double J;

  //Encode the non-zero matrix elements of the hamiltonian
  //on a given state using flipsh_ and mel_
  //i.e. all the state' such that <state'|H|state> = mel_(state') \neq 0
  //state' is encoded as the sequence of spin flips (saved as flips_) to be performed on state
  std::vector<std::vector<int> > flipsh_;
  std::vector<std::complex<double> > mel_;

public:

  Tc(int l,double beta_z=0,double beta_x=0,double beta_y=0, double J_vp=1):L(l),h_z(beta_z),h_x(beta_x),h_y(beta_y),J(J_vp){
   Init();}

  void Init(){
    mel_.resize(L*L+1+2*L*L);
    flipsh_.resize(L*L+1+2*L*L);

    //vertex terms sum A_v
    for (int v=0; v<L*L; v++){
      mel_[v+1]=-1.0*J;
      get_indices_star(v);
      flipsh_[v+1]=std::vector<int>(4);
      flipsh_[v+1][0]=l;
      flipsh_[v+1][1]=o;
      flipsh_[v+1][2]=r;
      flipsh_[v+1][3]=u;
    }
  
    //fields in x-direction sum h_x sigma_x
    for (int i=0; i<2*L*L; i++){
      int index0=L*L+1+i;
      mel_[index0]=h_x;
      flipsh_[index0]=std::vector<int>(1);
      flipsh_[index0][0]=i;
    }
  }


  //returns index of kronecker product
  inline int get_index1D(int x,int y,int AB){
    int z;
    z=L*2*x+y*2+AB;
    return z;
  }


  //saves indices of spins around a vertex v using the member variables l,o,r,u
  inline void get_indices_star(int v){
    int x,y;
    x=int(v/(L));
    y=v-L*x;
    int ym1,xp1;
    ym1=y-1;
    xp1=x+1;
    if (y==0){
       ym1=L-1;}
    if (x==(L-1)){
       xp1=0;} 
    l=get_index1D(x,ym1,1);
    o=get_index1D(x,y,0);
    r=get_index1D(x,y,1);
    u=get_index1D(xp1,y,0);
  }

  //saves indices of spins around a vertex v into the variables l1D,o1D,r1D,u1D
  inline void get_indices_star(int v,int &l1D,int &o1D, int &r1D,int &u1D){
    int x,y;
    x=int(v/(L));
    y=v-L*x;
    int ym1,xp1;
    ym1=y-1;
    xp1=x+1;
    if (y==0){
       ym1=L-1;}
    if (x==(L-1)){
       xp1=0;}
    l1D=get_index1D(x,ym1,1);
    o1D=get_index1D(x,y,0);
    r1D=get_index1D(x,y,1);
    u1D=get_index1D(xp1,y,0);
  }


  //saves indices of spins around a plaquette p using the member variables l,o,r,u
  inline void get_indices_plaquette(int p){
    int x=int(p)/int(L);
    int y=p-L*x;
    int yp1=y+1;
    int xm1=x-1;
    if (y==(L-1))
       	yp1=0;
    if (x==0)
       xm1=L-1;
    l=get_index1D(x,y,0);
    o=get_index1D(xm1,y,1);
    r=get_index1D(x,yp1,0);
    u=get_index1D(x,y,1);
  }


  //Finds the non-zero matrix elements of the hamiltonian
  //on the given state
  //i.e. all the state' such that <state'|H|state> = mel(state') \neq 0
  //state' is encoded as the sequence of spin flips to be performed on state
  void FindConn(const std::vector<int> & state, std::vector<std::vector<int> > & flipsh,std::vector<std::complex<double> > & mel){
    
    mel_[0]=0;

    //plaquette terms sum B_p
    for (int p=0;p<L*L; p++){
      get_indices_plaquette(p);
      std::complex<double> plaquette_cond=state[l]*state[o]*state[r]*state[u];
      mel_[0]= mel_[0]-plaquette_cond*J;
    }

    //fields in z-direction sum h_z sigma_z
    for (int i=0; i<2*L*L; i++){
      mel_[0]=mel_[0]+h_z*state[i];}

    //fields in y-direction sum h_y sigma_y
    for (int i=0; i<2*L*L; i++){
      int index0=L*L+1+i;
      mel_[index0].imag(h_y*state[i]);
    }   
       
    mel=mel_;
    flipsh=flipsh_;
  }

  int getL()const{
    return L;
  }

};


