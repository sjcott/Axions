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

const double tol = 1e-6;

const double eps = 0.1;
const double lambda = dz*(nz-1)/(1-0.25*eps*eps); // This essentially wavelength of sine string (not related to model parameters)

const double pi = 4*atan(1);
const double gr = 0.5*(1+sqrt(5));

int main(){

	Array phi1(2,nx,ny,nz,0.0), phi2(2,nx,ny,nz,0.0), A(3,nx,ny,nz,0.0), SOR_Fields(SORnx,3,0.0);

	double sigma[nx][nz], Fsigma[2], x, y, z, tolTest, Omega, distance, phi1Mag, phi2Mag, AMag, normal_dist, xs, zs, xs_sigma, zs_sigma, xs_sigma2, zs_sigma2,
		   paraVecMag, paraVecMag_sigma, a, b, sigma1, sigma2, distanceSqr1, distanceSqr2,finaldist;
	int i,j,k, pClosest;

	struct timeval start, end;
    gettimeofday(&start, NULL);

    string file_path = __FILE__;
    string dir_path = file_path.substr(0,file_path.find_last_of('/'));

    string SOR_inputPath = dir_path + "/Data/SOR_MooreFields.txt";
    string icPath = dir_path + "/Data/moore_ic.txt";
    //string testPath = dir_path + "/test.txt";

    ifstream SOR_input (SOR_inputPath.c_str());
    ofstream ic (icPath.c_str());
    //ofstream test (testPath.c_str());

    for(i=0;i<SORnx;i++){

    	SOR_input >> SOR_Fields(i,0) >> SOR_Fields(i,1) >> SOR_Fields(i,2);

    }

    // Define epsilon and Omega based off amplitude and wavelength. Make sure that epsilon is small for accuracy. Maybe add check that it's at least less than 1.

    Omega = 2*pi/lambda;

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

    for(k=z0;k<nz;k++){

    	z = (k-z0)*dz;

    	sigmaInit = z;

    	for(i=0;i<nx;i++){

    		x = (i-x0)*dx;

    		//attempt at 1st order initial guess approximation

    		//sigmaInit = z + eps*x*sin(Omega*z)/(1 - Omega*cos(Omega*z));

    		//sigma[i][k] = sigmaInit;

    		//cout << x << " " << z << endl;

    		tolTest=1;
    		int iterNum=0;

    		a = 0;
    		b = pi/Omega;

    		sigma1 = b - (b-a)/gr;
   			sigma2 = a + (b-a)/gr;

   			//cout << sigma1 << " " << sigma2 << " " << a << endl;

   			distanceSqr1 = pow(x - eps*cos(Omega*sigma1)/Omega,2) + pow(z - ( (1-0.25*eps*eps)*Omega*sigma1 + 0.125*eps*eps*sin(2*Omega*sigma1) )/Omega,2);
   			distanceSqr2 = pow(x - eps*cos(Omega*sigma2)/Omega,2) + pow(z - ( (1-0.25*eps*eps)*Omega*sigma2 + 0.125*eps*eps*sin(2*Omega*sigma2) )/Omega,2);

    		while(tolTest>tol){

    			/////////////////////// Newton - Rhapson method /////////////////////////////////////////////////////////

    			// tolTest = 0;

    			// // Minimises distance from string (ignoring y dimension)

   				// // Fsigma[0] = 2*eps*sin(sigma[i][k])*(x - eps*cos(sigma[i][k])/Omega) - 2*(1 - 0.25*eps*eps*(1 + cos(2*sigma[i][k])))*(z
   				// // 		  - ( (1-0.25*eps*eps)*sigma[i][k] - 0.125*eps*eps*sin(2*sigma[i][k]) )/Omega);

   				// // Fsigma[1] = 2*eps*cos(sigma[i][k])*(x - eps*cos(sigma[i][k])/Omega) + 2*pow(eps*sin(sigma[i][k]),2)/Omega
   				// //    		  - eps*eps*sin(2*sigma[i][k])*(z - ( (1-0.25*eps*eps)*sigma[i][k] - 0.125*eps*eps*sin(2*sigma[i][k]) )/Omega)
   				// //    		  + 2*pow(1-0.25*eps*eps*(1+cos(2*sigma[i][k])),2)/Omega;

   				// // String position, related quantities and derivatives wrt sigma

   				// xs = eps*cos(Omega*sigma[i][k])/Omega;
   				// zs = ( (1-0.25*eps*eps)*Omega*sigma[i][k] + 0.125*eps*eps*sin(2*Omega*sigma[i][k]) )/Omega;

   				// xs_sigma = -eps*sin(Omega*sigma[i][k]);
   				// zs_sigma = 1-0.25*eps*eps + 0.25*eps*eps*cos(2*Omega*sigma[i][k]);

   				// xs_sigma2 = -eps*Omega*cos(Omega*sigma[i][k]);
   				// zs_sigma2 = -0.5*eps*eps*Omega*sin(2*Omega*sigma[i][k]);

   				// paraVecMag = sqrt(pow(xs_sigma,2) + pow(zs_sigma,2));

   				// // paraVecMag_sigma = 0.5*( 2*xs_sigma*xs_sigma2 + 2*zs_sigma*zs_sigma2 )/paraVecMag;

   				// // Newton-Rhapson method for component of distance parallel to string

   				// // Fsigma[0] = ( xs_sigma*(x-xs) + zs_sigma*(z-zs) )/paraVecMag;

   				// // Fsigma[1] = ( (xs_sigma2*paraVecMag - xs_sigma*paraVecMag_sigma)*(x - xs) + (zs_sigma2*paraVecMag - zs_sigma*paraVecMag_sigma)*(z - zs) )/pow(paraVecMag,2)
   				// // 		  - ( pow(xs_sigma,2) + pow(zs_sigma,2) )/paraVecMag;

   				// // Newton-Rhapson method for distance (assumes that minimum distance is perpendicular to string)

   				// // Drives derivative of the distance squared to zero. Factor of 2 has been removed since it makes no difference to the zero point.

   				// Fsigma[0] = -xs_sigma*(x-xs) - zs_sigma*(z-zs);

   				// Fsigma[1] = -xs_sigma2*(x-xs) + pow(xs_sigma,2) - zs_sigma2*(z-zs) + pow(zs_sigma,2);


   				// sigma[i][k] += -Fsigma[0]/Fsigma[1];

   				// //tolTest = abs(Fsigma[0]/Fsigma[1]);
   				// if(abs(Fsigma[0])>tolTest){ tolTest = abs(Fsigma[0]); }

   				// iterNum += 1;

   				// if(iterNum > 1000000){

   				// 	cout << x << " " << z << " " << Fsigma[0] << " " << Fsigma[1] << " " << sigma[i][k] << endl;

   				// 	// iterNum = 0;

   				// 	// sigma[i][k] = rand

   				// }


   				/////////////////////////// Golden section search ///////////////////////////////////////////////

   				if(distanceSqr1<distanceSqr2){

   					b = sigma2;

   					//cout << "1>2 and " << a + gr*(b-a) << " " << sigma2 << endl;

   					sigma1 = b - (b-a)/gr;
   					sigma2 = a + (b-a)/gr;
   					
   					distanceSqr1 = pow(x - eps*cos(Omega*sigma1)/Omega,2) + pow(z - ( (1-0.25*eps*eps)*Omega*sigma1 + 0.125*eps*eps*sin(2*Omega*sigma1) )/Omega,2);
   					distanceSqr2 = pow(x - eps*cos(Omega*sigma2)/Omega,2) + pow(z - ( (1-0.25*eps*eps)*Omega*sigma2 + 0.125*eps*eps*sin(2*Omega*sigma2) )/Omega,2);

   				} else{

   					a = sigma1;

   					//cout << "2>1 and " << b - gr*(b-a) << " " << sigma1 << endl;

   					sigma1 = b - (b-a)/gr;
   					sigma2 = a + (b-a)/gr;

   					distanceSqr1 = pow(x - eps*cos(Omega*sigma1)/Omega,2) + pow(z - ( (1-0.25*eps*eps)*Omega*sigma1 + 0.125*eps*eps*sin(2*Omega*sigma1) )/Omega,2);
   					distanceSqr2 = pow(x - eps*cos(Omega*sigma2)/Omega,2) + pow(z - ( (1-0.25*eps*eps)*Omega*sigma2 + 0.125*eps*eps*sin(2*Omega*sigma2) )/Omega,2);

   				}

   				tolTest = b-a;

   				if(tolTest<=tol){

   					if(distanceSqr1>distanceSqr2){

   						sigma[i][k] = sigma2;

   						finaldist = distanceSqr2;


   					} else{

   						sigma[i][k] = sigma1;

   						finaldist = distanceSqr1;

   					}

   				}

    		}

    		// Use this value of sigma as the next initial guess

    		//sigmaInit = sigma[i][k];

    		// Test if the distance is equal to the dot product with unit normal

    		xs = eps*cos(Omega*sigma[i][k])/Omega;
   			zs = ( (1-0.25*eps*eps)*Omega*sigma[i][k] + 0.125*eps*eps*sin(2*Omega*sigma[i][k]) )/Omega;

   			xs_sigma = -eps*sin(Omega*sigma[i][k]);
   			zs_sigma = 1-0.25*eps*eps + 0.25*eps*eps*cos(2*Omega*sigma[i][k]);

			paraVecMag = sqrt(pow(xs_sigma,2) + pow(zs_sigma,2));

    		normal_dist = ( zs_sigma*(x-xs) - xs_sigma*(z-zs) )/paraVecMag;

     		//test << x << " " << z << " " << normal_dist << " " << sqrt(pow(x-xs,2) + pow(z-zs,2)) << " " << sqrt(finaldist) << endl;


    		for(j=0;j<ny;j++){

    			y = (j-y0)*dy;

    			distance = sqrt(pow(x-xs,2) + pow(y,2) + pow(z-zs,2));

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


     			if(normal_dist==0){

     				// To prevent division by zero

     				phi1(0,i,j,k) = 0;
     				phi2(0,i,j,k) = 0;
     				A(1,i,j,k) = 0;

     			} else{

     				phi1(0,i,j,k) = phi1Mag*normal_dist/distance;
     				phi2(0,i,j,k) = phi2Mag*normal_dist/distance;
     				A(1,i,j,k) = AMag*normal_dist/pow(distance,2);

     			}

     			if(y==0){

     				// To prevent division by zero

     				phi1(1,i,j,k) = 0;
     				phi2(1,i,j,k) = 0;
     				A(0,i,j,k) = 0;
     				A(2,i,j,k) = 0;

     			} else{

     				phi1(1,i,j,k) = phi1Mag*y/distance;
     				phi2(1,i,j,k) = phi2Mag*y/distance;
     				A(0,i,j,k) = -AMag*zs_sigma*y/(pow(distance,2)*paraVecMag);
     				A(2,i,j,k) = AMag*xs_sigma*y/(pow(distance,2)*paraVecMag);

     			}

     			// As long as this is not the first point in z (due to periodic boundary conditions) assign to z>0 using reflection symmetry across z=0

     			// if(k!=0){

     			// 	phi(0,i,j,nz-k) = phi(0,i,j,k);
     			// 	phi(1,i,j,nz-k) = phi(1,i,j,k);

     			// }

     			// Reflection across z=0

     			phi1(0,i,j,nz-1-k) = phi1(0,i,j,k);
     			phi1(1,i,j,nz-1-k) = phi1(1,i,j,k);

     			phi2(0,i,j,nz-1-k) = phi2(0,i,j,k);
     			phi2(1,i,j,nz-1-k) = phi2(1,i,j,k);

     			A(0,i,j,nz-1-k) = A(0,i,j,k);
     			A(1,i,j,nz-1-k) = A(1,i,j,k);
     			if(k!=nz-1){ A(2,i,j,nz-2-k) = -A(2,i,j,k); }

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


