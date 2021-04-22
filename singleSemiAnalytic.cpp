#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <sys/time.h>
#include <omp.h>
#include <complex>

#include "array.hpp"

using namespace std;

// Calculated the radial profile of global or Abelian-Higgs, straight string for radiation study. Both parameters of the global model can be removed by rescalings so here we are solving
// a system with no free parameters. |Phi| is proportional to f_a and the length scale x is proportional to lambda^{-1/2}. Here we use the rescaled phi which has f_a=1 and lambda=2.

// Abelian-Higgs introduces the gauge coupling parameter.


const int nx =30001;
const double h=0.01;
const double w=1.5; // Relaxation factor > 1. Trial and error to get better convergence.
const int maxIter = 2000000;
const double tol = 1e-3;
const int ic_type = 0;


const double pi = 4*atan(1);

const double n = 1; // Winding of the string
const double lambda = 1;
const double eta = 1;
const double g = 0; // Gauge coupling (g=1 is the BPS limit). Set to zero to return to a global string

int main(){

	struct timeval start, end;
    gettimeofday(&start, NULL);

    string file_path = __FILE__;
    string dir_path = file_path.substr(0,file_path.find_last_of('/'));

    string SOR_FieldsPath = dir_path + "/Data/SOR_Fields.txt";

    ifstream ic (SOR_FieldsPath.c_str());

    // Initialise the field array

	double tolTest;

	int i,j,iterNum;

	Array phi(nx,0.0), Fphi(nx,2,0.0), A(nx,0.0), FA(nx,2,0.0);

	if(ic_type==0){

		for(i=0;i<nx;i++){

			phi(i) = eta*tanh(0.5*h*i);
			A(i)  = 0;

		}

	} else if(ic_type==1){

		for(i=0;i<nx;i++){

			ic >> phi(i) >> A(i);

		}

	}

    ofstream SOR_Fields (SOR_FieldsPath.c_str());

	// Set boundary conditions
	phi(0) = 0;

	A(0) = 0;

	tolTest = 1;
	iterNum = 0;
	while(tolTest>tol and iterNum<maxIter){

		tolTest = 0;
		iterNum += 1;

		if((100*iterNum)%maxIter==0){

			cout << "\rIteration number: " << iterNum << flush;

		}

		if(iterNum==maxIter){

			cout << "\rMaximum number of iterations reached" << flush;

		}

		// Set derivative boundary conditions

		phi(nx-1) = phi(nx-2);
		A(nx-1) = A(nx-2);


		// Looping over all positions in grid except boundaries.
		# pragma omp parallel for default(none) shared(phi,A,Fphi,FA)
		for(j=1;j<nx-1;j++){


			Fphi(j,0) = ( pow((n-g*A(j))/(j*h),2) + 0.5*lambda*(pow(phi(j),2) - pow(eta,2)) )*phi(j) - ( phi(j+1) - phi(j-1) )/(2*j*h*h) - ( phi(j+1) - 2*phi(j) + phi(j-1) )/(h*h);

			Fphi(j,1) = pow((n-g*A(j))/(j*h),2) + 0.5*lambda*(3*pow(phi(j),2) - pow(eta,2)) + 2/(h*h);  // Derivative of Fphi[0] with respect to phi[i]


			FA(j,0) = -2*g*(n-g*A(j))*pow(phi(j),2) + ( A(j+1) - A(j-1) )/(2*j*h*h) - ( A(j+1) - 2*A(j) + A(j-1) )/(h*h);

			FA(j,1) = 2*pow(g*phi(j),2) + 2/(h*h);


			// Now update the field values by SOR.

			phi(j) += - w*Fphi(j,0)/Fphi(j,1);
			A(j) += -w*FA(j,0)/FA(j,1);

		}

		for(j=1;j<nx-1;j++){

			if(abs(Fphi(j,0))>tolTest){ tolTest = abs(Fphi(j,0)); }
			if(abs(FA(j,0))>tolTest){ tolTest = abs(FA(j,0)); }

		}

	}

	for(j=0;j<nx;j++){

		SOR_Fields << phi(j) << " " << A(j) << endl;

	}

	cout << "\rNumber of iterations: " << iterNum << endl;

	gettimeofday(&end, NULL);

    cout << "Time taken: " << end.tv_sec - start.tv_sec << "s" << endl;


	return 0;

}