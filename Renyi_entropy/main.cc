//Author: Agnes Valenti. 

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
#include <complex>
#include "sampler_renyi.cc"
#include <sys/time.h>

int main(){

    //lattice size (N=2*L*L spins)
    int L=3;

    //n_hvalues-start_hi=number of field strength steps
    int n_hvalues=60;
    int start_h_i=0;

    //field in (hx,hy,hz)-direction
    double hx,hy,hz;

    //choose the field direction:
    //  field_direction=1 for the self-dual line h(1,0,1)
    //  field_direction=2 for h(1,0.2,0.5)
    //  field_direction=3 for (h,0.2,h)
    int field_direction=1;

    //for ground states, exc=false. To compute
    //energies of the first excited state choose exc=true.
    //(obtained via orthogonalization added to the cost function)
    //only compatible with field_direction=2
    bool exc=false;
    if (exc){
      if(field_direction!=2){
        std::cout<<"ERROR: Invalid combination of field direction and choice of excited state."<<std::endl;
        std::cout<<"For the computation of excited state energies, choose field direction 2."<<std::endl;
          std::cout<<std::endl;
          std::abort();
      }
      start_h_i=10;
    }
    if (field_direction==1){
      hx=1;
      hy=0;
      hz=1;
    }
    if (field_direction==2){
      hx=1;
      hy=0.2;
      hz=0.5;
    }

    if (field_direction==3){
      hx=1;
      hy=0.2;
      hz=1;
    }



    std::cout<<"Number of spins: "<<2*L*L<<std::endl;
    std::cout<<std::endl;

    std::cout<<"Number of Monte Carlo sweeps: "<<nsweeps<<std::endl;
    std::cout<<std::endl;

    if (exc){
      std::cout<<"Computing the cRBM second Renyi entropies S_2 of the first excited state for the field h in direction ("<<hx<<", "<<hy<<", "<<hz<<")"<<std::endl;
        std::cout<<"in the range 0.1<h<0.61 in steps of size delta_h=0.01..."<<std::endl;
        std::cout<<std::endl;
    }
    else{
        std::cout<<"Computing the cRBM second Renyi entropies S_2 of the ground state for the field h in direction ("<<hx<<", "<<hy<<", "<<hz<<")"<<std::endl;
        std::cout<<"in the range 0<h<0.61 in steps of size delta_h=0.01..."<<std::endl;
        std::cout<<std::endl;
    }

    //save entropies into file
    std::ofstream file_b;
    std::ostringstream fileNameStream_b("");
    fileNameStream_b<<"entropies.txt";
    std::string fileName_b=fileNameStream_b.str();
    file_b.open(fileName_b.c_str());


    //uncomment in order to parallelize
    //#pragma omp parallel for num_threads(n_hvalues-start_h_i)
    for (int h_i=start_h_i; h_i<n_hvalues; h_i++){
     
     //field strength
     double h=0.01+h_i*0.01; //h_i*0.1-50.0;
     
     //load translational symmetries
     std::ostringstream fileNameStream_symm("");
     fileNameStream_symm<<"symmetries_L"<<L<<".txt";
     std::string fileName_symm=fileNameStream_symm.str();

       //load cRBM weights
      std::ostringstream fileNameStream_weights("");
      if (exc){
         fileNameStream_weights<<"../weights/L3/exc_field_hx"<<hx<<"_hy"<<hy<<"_hz"<<hz<<"/weights"<<h_i<<".txt";

      }
      else{
         fileNameStream_weights<<"../weights/L3/field_hx"<<hx<<"_hy"<<hy<<"_hz"<<hz<<"/weights"<<h_i<<".txt";
      }
      std::string fileName_weights=fileNameStream_weights.str();
  
      
      //create cRBM objects (2 copies) to sample from
      cRBM Wf1(fileName_weights,fileName_symm); //----------------------------L!
      cRBM Wf2(fileName_weights,fileName_symm); //----------------------------L!

        
      Sampler_Renyi<cRBM> sampler1(Wf1, Wf2,L,-1);
    
      
      //run MC sampling to calculate energies 
      sampler1.Run_MC();
 

      //entropy, var=MC sampling error
      double entropy, var;
      sampler1.OutputEntropy(entropy,var);
      entropy=1.0/(1.0-2.0)*std::log2(entropy);

      file_b<<h<<" "<<entropy<<std::endl;
      
      std::cout<<"Field strength: h="<<h<<". Computed cRBM second Renyi entropy S_2="<<entropy<<" +/- "<<var<<std::endl;  //put here MC error
      std::cout<<std::endl;

    }
    
        file_b.close();


    return 0;
    }

 

