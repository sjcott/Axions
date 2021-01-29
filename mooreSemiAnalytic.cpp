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


const int nx =20001;
const double h=0.01;
const double w=1.5; // Relaxation factor > 1. Trial and error to get better convergence.
const int maxIter = 2000000;
const double tol = 1e-6;
const int ic_type = 0;


const double pi = 4*atan(1);

const double n = 1; // Winding of string 1
const double m = 1; // Winding of string 2

// Model parameters

const double eta1 = 1; // Can fix without loss of generality (rescales solutions)
const double eta2 = 1;
const double lambda1 = 2; // Can fix without loss of generality (rescales solutions)
const double lambda2 = 2;
const double lambda12 = 0;
const int q1 = 2;
const int q2 = q1-1; // Moore model uses q2 = q1-1 with q1 integer so that (for (1,1) strings) the global charge of the string is 1.
const double g = sqrt(0.5*lambda1*pow(eta1,2)/(q1*q1 + q2*q2)); // Specific choice that sets all masses equal (assuming m1=m2)

const double rmin = pi/2; // Lower cut off radius. Here chosen (as by Moore) to be pi/mass

int main(){

	struct timeval start, end;
    gettimeofday(&start, NULL);

    string file_path = __FILE__;
    string dir_path = file_path.substr(0,file_path.find_last_of('/'));

    string SOR_MooreFieldsPath = dir_path + "/Data/SOR_MooreFields.txt";
    string propsPath = dir_path + "/Data/props.txt";

    ifstream ic (SOR_MooreFieldsPath.c_str());

    // Initialise the field array

	double tolTest;

	int i,j,iterNum;

	Array phi1(nx,0.0), Fphi1(nx,2,0.0), phi2(nx,0.0), Fphi2(nx,2,0.0), A(nx,0.0), FA(nx,2,0.0);

	if(ic_type==0){

		for(i=0;i<nx;i++){

			phi1(i) = tanh(0.5*h*i);
			phi2(i) = tanh(0.5*h*i);
			A(i)  = 0;

		}

	} else if(ic_type==1){

		for(i=0;i<nx;i++){

			ic >> phi1(i) >> phi2(i) >> A(i);

		}

	}

    ofstream SOR_Fields (SOR_MooreFieldsPath.c_str());
    ofstream props (propsPath.c_str());

	// Set boundary conditions
	phi1(0) = 0;
	phi2(0) = 0;
	A(0) = 0;

	// phi1(nx-1) = eta1;
	// phi2(nx-1) = eta2;
	// A(nx-1) = (n*q1*pow(eta1,2) + m*q2*pow(eta2,2))/(g*(pow(q1*eta1,2)+pow(q2*eta2,2)));

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

		phi1(nx-1) = phi1(nx-2);
		phi2(nx-1) = phi2(nx-2);
		A(nx-1) = A(nx-2);


		// Looping over all positions in grid except boundaries.
		# pragma omp parallel for default(none) shared(phi1,phi2,A,Fphi1,Fphi2,FA)
		for(j=1;j<nx-1;j++){


			Fphi1(j,0) = ( pow( (n-q1*g*A(j))/(j*h) ,2) + 0.5*lambda1*( pow(phi1(j),2) - pow(eta1,2) ) + 0.5*lambda12*( pow(phi2(j),2) - pow(eta2,2) ) )*phi1(j)
					   - ( phi1(j+1) - phi1(j-1) )/(2*j*h*h) - ( phi1(j+1) - 2*phi1(j) + phi1(j-1) )/(h*h);

			Fphi1(j,1) = pow( (n-q1*g*A(j))/(j*h) ,2) + 0.5*lambda1*( 3*pow(phi1(j),2) - pow(eta1,2) ) + 0.5*lambda12*( pow(phi2(j),2) - pow(eta2,2) ) + 2/(h*h);


			Fphi2(j,0) = ( pow( (m-q2*g*A(j))/(j*h) ,2) + 0.5*lambda2*( pow(phi2(j),2) - pow(eta2,2) ) + 0.5*lambda12*( pow(phi1(j),2) - pow(eta1,2) ) )*phi2(j)
					   - ( phi2(j+1) - phi2(j-1) )/(2*j*h*h) - ( phi2(j+1) - 2*phi2(j) + phi2(j-1) )/(h*h);

			Fphi2(j,1) = pow( (m-q2*g*A(j))/(j*h) ,2) + 0.5*lambda2*( 3*pow(phi2(j),2) - pow(eta2,2) ) + 0.5*lambda12*( pow(phi1(j),2) - pow(eta1,1) ) +2/(h*h);


			FA(j,0) = - 2*q1*g*(n-q1*g*A(j))*pow(phi1(j),2) - 2*q2*g*(m-q2*g*A(j))*pow(phi2(j),2) + ( A(j+1) - A(j-1) )/(2*j*h*h) - ( A(j+1) - 2*A(j) + A(j-1) )/(h*h);

			FA(j,1) = 2*pow(q1*g*phi1(j),2) + 2*pow(q2*g*phi2(j),2) + 2/(h*h);


			// Now update the field values by SOR.

			phi1(j) += - w*Fphi1(j,0)/Fphi1(j,1);
			phi2(j) += - w*Fphi2(j,0)/Fphi2(j,1);
			A(j) += -w*FA(j,0)/FA(j,1);

		}

		for(j=1;j<nx-1;j++){

			if(abs(Fphi1(j,0))>tolTest){ tolTest = abs(Fphi1(j,0)); }
			if(abs(Fphi2(j,0))>tolTest){ tolTest = abs(Fphi2(j,0)); }
			if(abs(FA(j,0))>tolTest){ tolTest = abs(FA(j,0)); }

		}

	}

	for(j=0;j<nx;j++){

		SOR_Fields << phi1(j) << " " << phi2(j) << " " << A(j) << endl;

	}

	// Calculate the energy and effective kappa

	double E = 0;

	for(j=1;j<nx-1;j++){

		E += 2*pi*h*j*h*( pow((phi1(j+1)-phi1(j-1))/(2*h),2) + pow((phi2(j+1)-phi2(j-1))/(2*h),2) + pow((n-q1*g*A(j))*phi1(j)/(j*h),2) + pow((m-q2*g*A(j))*phi2(j)/(j*h),2)
		   + 0.5*pow((A(j+1)-A(j-1))/(2*j*h*h),2) + 0.25*lambda1*pow(pow(phi1(j),2)-pow(eta1,2),2) + 0.25*lambda2*pow(pow(phi2(j),2)-pow(eta2,2),2)
		   + 0.5*lambda12*(pow(phi1(j),2)-pow(eta1,2))*(pow(phi2(j),2)-pow(eta2,2)) );

	}

	double fa = 2*pow(eta1*eta2,2)/(pow(q1*eta1,2) + pow(q2*eta2,2)); // This is goldstone mode decay constant squared

	double kappa = E/(pi*fa) - log((nx-1)*h/rmin);	

	props << E << " " << kappa << endl;

	cout << "\rNumber of iterations: " << iterNum << endl;

	gettimeofday(&end, NULL);

    cout << "Time taken: " << end.tv_sec - start.tv_sec << "s" << endl;


	return 0;

}