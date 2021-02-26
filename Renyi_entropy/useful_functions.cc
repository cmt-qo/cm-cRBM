//Author: Agnes Valenti
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <string>
#include <sstream>
#include <random>
#include <iomanip>
#include <limits>
#include <ctime>


inline int get_index1D_(int x,int y,int AB, int L){
    //returns index of kronecker product
    int z;
    z=L*2*x+y*2+AB;
    return z;
    }

//returns (as call by reference variables s1, s2) the indices of the two vertices adjacent to spin
inline void get_2vertices(int spin, int & s1, int & s2, int L_){
   int spin_x,spin_y,spin_AB;
   spin_x=spin/(2*L_);
   spin_y=(spin-spin_x*2*L_)/2;
   spin_AB=(spin-spin_x*2*L_-spin_y*2);
   int xs1,xs2,ys1,ys2;
   if (spin_AB==0){
      ys1=spin_y;
      ys2=spin_y;
      if (spin_x==0)
        xs1=L_-1;
      else
        xs1=(spin_x-1);
      xs2=spin_x;
      }
   if (spin_AB==1){
      xs1=spin_x;
      xs2=spin_x;
      ys1=spin_y;
      ys2=(spin_y+1)%L_;
      }
   s1=L_*xs1+ys1;
   s2=L_*xs2+ys2;
   }

//returns (as call by reference variables p1, p2) the indices of the two plaquettes adjacent to spin
inline void get_2plaquettes(int spin, int & p1, int & p2, int L_){
   int spin_x,spin_y,spin_AB;
   spin_x=spin/(2*L_);
   spin_y=(spin-spin_x*2*L_)/2;
   spin_AB=(spin-spin_x*2*L_-spin_y*2);
   int xp1,xp2,yp1,yp2;
   if (spin_AB==1){
      yp1=spin_y;
      yp2=spin_y;
      xp1=spin_x;
      xp2=(spin_x+1)%L_;
      }
   if (spin_AB==0){
      xp1=spin_x;
      xp2=spin_x;
      yp1=spin_y;
      if (spin_y==0)
        yp1=L_-1;
      else
        yp1=(spin_y-1);
      yp2=spin_y%L_;
      }
   p1=L_*xp1+yp1;
   p2=L_*xp2+yp2;
   }


inline void get_indices_star_(int v,int &l1D,int &o1D, int &r1D,int &u1D, int L){
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
    l1D=get_index1D_(x,ym1,1,L);
    o1D=get_index1D_(x,y,0,L);
    r1D=get_index1D_(x,y,1,L);
    u1D=get_index1D_(xp1,y,0,L);
  }

inline void get_indices_plaquette_(int p,int &l1D,int &o1D, int &r1D,int &u1D, int L){
    int x=int(p)/int(L);
    int y=p-L*x;
    int yp1=y+1;
    int xm1=x-1;
    if (y==(L-1))
       	yp1=0;
    if (x==0)
       xm1=L-1;
    l1D=get_index1D_(x,y,0,L);
    o1D=get_index1D_(xm1,y,1,L);
    r1D=get_index1D_(x,yp1,0,L);
    u1D=get_index1D_(x,y,1,L);
    }

