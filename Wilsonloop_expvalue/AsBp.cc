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

//constructs W(m=4) operator
class AsBp{

  //lattice size
  const int L;
 
  //useful variables to save spin indices
  int l,o,r,u;

  double J;
  
  //Encode the non-zero matrix elements of the operator W(m=4)
  //on a given state using flipsh_ and mel_
  //i.e. all the state' such that <state'|W(m=4)|state> = mel_(state') \neq 0
  //state' is encoded as the sequence of spin flips (saved as flips_) to be performed on state
  std::vector<std::vector<int> > flipsh_;
  std::vector<std::complex<double> > mel_;
public:

  AsBp(int l,double J_vp=1):L(l),J(J_vp){
    Init();}

  void Init(){
    mel_.resize(1);
    flipsh_.resize(1);
    flipsh_[0]=std::vector<int>(8);
    flipsh_[0][0]=0;
    flipsh_[0][1]=2;
    flipsh_[0][2]=3;
    flipsh_[0][3]=2*L-1;
    flipsh_[0][4]=2*L+3;
    flipsh_[0][5]=4*L-1;
    flipsh_[0][6]=4*L;
    flipsh_[0][7]=4*L+2;
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

  //Finds the non-zero matrix elements of the operator W(m=4)
  //on the given state
  //i.e. all the state' such that <state'|W(m=4)|state> = mel(state') \neq 0
  //state' is encoded as the sequence of spin flips to be performed on state
  void FindConn(const std::vector<int> & state, std::vector<std::vector<int> > & flipsh,std::vector<std::complex<double> > & mel){
    std::vector<int> indices(4);
    indices[0]=0; 
    indices[1]=1; 
    indices[2]=L+1; 
    indices[3]=L+2; 
    mel_[0]=0;
    std::complex<double> plaquette_cond=1;
    for (int p=0;p<4; p++){
      get_indices_plaquette(indices[p]);
      plaquette_cond=plaquette_cond*double(state[l])*double(state[o])*double(state[r])*double(state[u]);
    }
    mel_[0]= mel_[0]-plaquette_cond*J;
       
       
    mel=mel_;
    flipsh=flipsh_;
  }

  int getL()const{
    return L;
  }

};


