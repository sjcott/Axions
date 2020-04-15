#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <sys/time.h>
#include <omp.h>
#include <complex>

#include "header.hpp"

using namespace std;

// Calculated the radial profile of a global, straight string for radiation study. Both parameters of the model can be removed by rescalings so here we are solving
// a system with no free parameters. |Phi| is proportional to f_a and the length scale x is proportional to lambda^{-1/2}. Here we use the rescaled phi.


const int nx =1000001;
const double h=0.001;
const double w=1.5; // Relaxation factor > 1. Trial and error to get better convergence.
const int maxIter = 20000000;
const double tol = 1e-3;
const int ic_type = 0;


const double pi = 4*atan(1);
const double n = 1; // Winding of the string

int main(){

	struct timeval start, end;
    gettimeofday(&start, NULL);

    string file_path = __FILE__;
    string dir_path = file_path.substr(0,file_path.find_last_of('/'));

    string SOR_FieldsPath = dir_path + "/SOR_Fields.txt";

    ifstream ic (SOR_FieldsPath.c_str());

    // Initialise the field array

	double tolTest;

	int i,j,iterNum;

	Array phi(nx,0.0), Fphi(nx,2,0.0);

	if(ic_type==0){

		for(i=0;i<nx;i++){

			phi(i) = tanh(0.5*h*i);

		}

	} else if(ic_type==1){

		for(i=0;i<nx;i++){

			ic >> phi(i);

		}

	}

    ofstream SOR_Fields (SOR_FieldsPath.c_str());

	// Set boundary conditions
	phi(0) = 0;
	phi(nx-1) = 1;

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


		// Looping over all positions in grid except boundaries.
		# pragma omp parallel for default(none) shared(phi,Fphi)
		for(j=1;j<nx-1;j++){


			Fphi(j,0) = ( pow(n/(j*h),2) + phi(j)*phi(j) - 1 )*phi(j) - ( phi(j+1) - phi(j-1) )/(2*j*h*h) - ( phi(j+1) - 2*phi(j) + phi(j-1) )/(h*h);

			Fphi(j,1) = pow(n/(j*h),2) + 3*phi(j)*phi(j) - 1 + 2/(h*h);  // Derivative of Fphi[0] with respect to phi[i]


			// Now update the field values by SOR.

			phi(j) = phi(j) - w*Fphi(j,0)/Fphi(j,1);

		}

		for(j=1;j<nx-1;j++){

			if(abs(Fphi(j,0))>tolTest){ tolTest = abs(Fphi(j,0)); }

		}

	}

	for(j=0;j<nx;j++){

		SOR_Fields << phi(j) << endl;

	}

	cout << "\rNumber of iterations: " << iterNum << endl;

	gettimeofday(&end, NULL);

    cout << "Time taken: " << end.tv_sec - start.tv_sec << "s" << endl;


	return 0;

}