#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <sys/time.h>
#include <omp.h>
#include <complex>

#include "array.hpp"

using namespace std;

// Perturb a straight, global string solution, evolve the system and investigate how it radiates.

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                           Class Definitions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// class Array{

//     public:

//         Array(int,int,int,int,double); // Constructor for 4d arrays
//         Array(int,int,int,double); // Constructor for 3d arrays

//         int A, B, C, D; // Sizes of each dimension
//         vector<double> array; //Empty vector container for the array

//         double& operator()(int a, int b, int c, int d){ // Overloads operator to allow easy indexing of 4d arrays

//             return array[D*(C*(B*a + b) + c) + d];

//         }

//         double& operator()(int a, int b, int c){ // Same as above for 3d arrays

//             return array[C*(B*a + b) + c];

//         }

// };

// // Defining array constructors

// Array::Array(int a, int b, int c, int d, double v){

//     A = a;
//     B = b;
//     C = c;
//     D = d;

//     array.assign(A*B*C*D, v);

// }

// Array::Array(int a, int b, int c, double v){

//     A = a;
//     B = b;
//     C = c;

//     array.assign(A*B*C,v);

// }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                 		  Parameters & Declarations
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

const int nx = 201;
const int ny = 201;
const int nz = 51;
const double dx = 0.7;
const double dy = 0.7;
const double dz = 0.7;

const int n = 1;
const int m = 1; // m and n are useful for now, code assumes they are both 1.
const int q1 = 2;
const int q2 = q1 - 1;

const double eta1 = 1;
const double eta2 = 1;
const double lambda1 = 2;
const double g = sqrt(0.5*lambda1*pow(eta1,2)/(q1*q1 + q2*q2)); // Specific choice that sets all masses equal (assuming m1=m2)

const int SORnx = 20001;
const double SORa = 0.01;

const double eps = 0;
const double lambda = dz*(nz-1)/(1-0.25*eps*eps); // This essentially wavelength of sine string (not related to model parameters)

const double pi = 4*atan(1);
const double gr = 0.5*(1+sqrt(5));

int main(){

	Array phi1(2,nx,ny,nz,0.0), phi2(2,nx,ny,nz,0.0), A(3,nx,ny,nz,0.0), SOR_Fields(SORnx,3,0.0);

	double x, y, distance, phi1Mag, phi2Mag, AMag;
	int i,j,k, pClosest;

	struct timeval start, end;
    gettimeofday(&start, NULL);

    string file_path = __FILE__;
    string dir_path = file_path.substr(0,file_path.find_last_of('/'));

    string SOR_inputPath = dir_path + "/SOR_MooreFields.txt";
    string icPath = dir_path + "/moore_ic.txt";
    //string testPath = dir_path + "/test.txt";

    ifstream SOR_input (SOR_inputPath.c_str());
    ofstream ic (icPath.c_str());
    //ofstream test (testPath.c_str());

    for(i=0;i<SORnx;i++){

    	SOR_input >> SOR_Fields(i,0) >> SOR_Fields(i,1) >> SOR_Fields(i,2);

    }

    // Grid positions of zero. Difference due to the fact that the z direction has periodic boundary conditions emposed. x=0 and y=0 are not on the grid.

    // double x0 = (nx+1)/2;
    // double y0 = (ny+1)/2;
    // int z0 = round(nz/2);

    int x0 = round((nx-1)/2);
    int y0 = round((ny-1)/2);
    int z0 = round((nz-1)/2);

    // Now need to perturb the string with ability to do it in multiple different ways.

    // Sine wave in x direction

    // Find minimum distance to string and treat this as radial direction for straight string solution.
    // Loop over z first so can use previous solution as initial guess for next.

    double sigmaInit = 0;

    for(k=0;k<nz;k++){
    	for(i=0;i<nx;i++){

    		x = (i-x0)*dx;

    		for(j=0;j<ny;j++){

    			y = (j-y0)*dy;

    			distance = sqrt(pow(x,2) + pow(y,2));

    			pClosest = round(distance/SORa);

    			if(pClosest==0){

    				// 1st order interpolation since only have grid points on one side.

     				phi1Mag = ( SOR_Fields(pClosest+1,0)*distance - SOR_Fields(pClosest,0)*(distance-SORa) )/SORa;
     				phi2Mag = ( SOR_Fields(pClosest+1,1)*distance - SOR_Fields(pClosest,1)*(distance-SORa) )/SORa;
     				AMag = ( SOR_Fields(pClosest+1,2)*distance - SOR_Fields(pClosest,2)*(distance-SORa) )/SORa;

     			} else if(pClosest<SORnx){

     				// 2nd order interpolation

     				phi1Mag = ( SOR_Fields(pClosest-1,0)*(distance-pClosest*SORa)*(distance-(pClosest+1)*SORa) - 2*SOR_Fields(pClosest,0)*(distance-(pClosest-1)*SORa)*(distance-(pClosest+1)*SORa) +
     						   SOR_Fields(pClosest+1,0)*(distance-(pClosest-1)*SORa)*(distance-pClosest*SORa) )/(2*SORa*SORa);

     				phi2Mag = ( SOR_Fields(pClosest-1,1)*(distance-pClosest*SORa)*(distance-(pClosest+1)*SORa) - 2*SOR_Fields(pClosest,1)*(distance-(pClosest-1)*SORa)*(distance-(pClosest+1)*SORa) +
     						   SOR_Fields(pClosest+1,1)*(distance-(pClosest-1)*SORa)*(distance-pClosest*SORa) )/(2*SORa*SORa);

     				AMag = ( SOR_Fields(pClosest-1,2)*(distance-pClosest*SORa)*(distance-(pClosest+1)*SORa) - 2*SOR_Fields(pClosest,2)*(distance-(pClosest-1)*SORa)*(distance-(pClosest+1)*SORa) +
     						 SOR_Fields(pClosest+1,2)*(distance-(pClosest-1)*SORa)*(distance-pClosest*SORa) )/(2*SORa*SORa);

     			} else{

     				phi1Mag = eta1;
     				phi2Mag = eta2;

     				if(g==0){ AMag = 0; }
     				else{ AMag = n/g; }

     				cout << "Off straight string solution grid" << endl;

     			}

     			// Now need to set phase of phi and split A_theta (theta defined with respect to string) into cartesian components


     			if(distance<0.1*min(dx,min(dy,dz))){

     				// To prevent division by zero

     				phi1(0,i,j,k) = 0;
     				phi2(0,i,j,k) = 0;
     				A(1,i,j,k) = 0;

     			} else{

     				phi1(0,i,j,k) = phi1Mag*x/distance;
     				phi2(0,i,j,k) = phi2Mag*x/distance;
     				A(1,i,j,k) = AMag*x/pow(distance,2);

     			}

     			if(distance<0.1*min(dx,min(dy,dz))){

     				// To prevent division by zero

     				phi1(1,i,j,k) = 0;
     				phi2(1,i,j,k) = 0;
     				A(0,i,j,k) = 0;
     				A(2,i,j,k) = 0;

     			} else{

     				phi1(1,i,j,k) = phi1Mag*y/distance;
     				phi2(1,i,j,k) = phi2Mag*y/distance;
     				A(0,i,j,k) = -AMag*y/pow(distance,2);
     				A(2,i,j,k) = 0;

     			}

     			// As long as this is not the first point in z (due to periodic boundary conditions) assign to z>0 using reflection symmetry across z=0

     			// if(k!=0){

     			// 	phi(0,i,j,nz-k) = phi(0,i,j,k);
     			// 	phi(1,i,j,nz-k) = phi(1,i,j,k);

     			// }

    		}

    	}
    }

    // Output field to file

    for(i=0;i<nx;i++){
    	for(j=0;j<ny;j++){
    		for(k=0;k<nz;k++){

    			// Convert gauge field to lattice link variable for use in evolution code

    			ic << phi1(0,i,j,k) << " " << phi1(1,i,j,k) << " " << phi2(0,i,j,k) << " " << phi2(1,i,j,k) << " " << dx*g*A(0,i,j,k) << " " << dy*g*A(1,i,j,k) << " " << dz*g*A(2,i,j,k) << endl;

    		}
    	}
    }

    // for(i=0;i<nx;i++){

    // 	x = (i-x0)*dx;

    // 	for(k=0;k<nz;k++){

    // 		z = (k-z0)*dz;
    // 		y = (j-y0)*dy;

    // 		for(j=0;j<nz;j++){

    // 			y = (j-y0)*dy;

    // 			distance = sqrt( pow(x-eps*cos(sigma[i][k])/Omega,2) + y*y + pow(z-( (1-0.25*eps*eps)*sigma[i][k] - 0.125*eps*eps*sin(2*sigma[i][k]) )/Omega,2) );

    // 			pClosest = round(distance/SORa);

    // 			if(pClosest==0){

    // 				// 1st order interpolation since only have grid points on one side.

    // 				phiMag = ( SOR_Field[pClosest+1]*distance - SOR_Field[pClosest]*(distance-SORa) )/SORa;

    // 			} else if(pClosest<SORnx){

    // 				// 2nd order interpolation

    // 				phiMag = ( SOR_Field[pClosest-1]*(distance-pClosest*SORa)*(distance-(pClosest+1)*SORa) - 2*SOR_Field[pClosest]*(distance-(pClosest-1)*SORa)*(distance-(pClosest+1)*SORa) +
    // 						   SOR_Field[pClosest+1]*(distance-(pClosest-1)*SORa)*(distance-pClosest*SORa) )/(2*SORa*SORa);

    // 			} else{

    // 				phiMag = 1;

    // 			}



    // 			// Now need to set phase of phi. Relevant angle is tan(alpha) = delta_y/sqrt(delta_x^2 + delta_z^2). Real part is phiMag*cos(alpha), imaginary is phiMag*sin(alpha)

    // 			// Test if found distance perpendicular to string with unit normal.

    // 			xs = eps*cos(sigma[i][k])/Omega;
   	// 			zs = ( (1-0.25*eps*eps)*sigma[i][k] - 0.125*eps*eps*sin(2*sigma[i][k]) )/Omega;

   	// 			xs_sigma = -eps*sin(sigma[i][k])/Omega;
   	// 			zs_sigma = ( (1-0.25*eps*eps) - 0.25*eps*eps*cos(2*sigma[i][k]) )/Omega;

    // 			normal_dist = ( (1/Omega)*( (1-0.25*eps*eps) - 0.25*eps*eps*cos(2*sigma[i][k]) )*(x - eps*cos(sigma[i][k])/Omega) + (eps/Omega)*sin(sigma[i][k])*( (1-0.25*eps*eps)*sigma[i][k] 
    // 						- 0.125*eps*eps*sin(2*sigma[i][k]) )/Omega )/sqrt( pow(eps*sin(sigma[i][k])/Omega,2) + pow(( (1-0.25*eps*eps) - 0.25*eps*eps*cos(2*sigma[i][k]) )/Omega,2) );

    // 			test << pow(normal_dist,2) + y*y << " " << distance << endl;

    // 			if(normal_dist==0){

    // 				Re_phi = 0;

    // 			} else{

    // 				Re_phi = phiMag*normal_dist/sqrt(pow(normal_dist,2)+y*y);

    // 			}

    // 			if(y==0){

    // 				Im_phi = 0;

    // 			} else{

    // 				Im_phi = -phiMag*y/sqrt(pow(normal_dist,2)+y*y);

    // 			}

    // 			ic << Re_phi << " " << Im_phi << endl;

    // 		}
    // 	}
    // }



    gettimeofday(&end,NULL);

    cout << "Time taken: " << end.tv_sec - start.tv_sec << "s" << endl;


	return 0;

}


