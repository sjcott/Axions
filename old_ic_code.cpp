#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <sys/time.h>
#include <omp.h>
#include <complex>
#include <random>

#include "array.hpp"

using namespace std;
typedef complex<double> dcmplx;

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

const string ic_type = "string"; // Which type of initial condition generation to use. ""string" bases initial conditions on straight string solutions (sine perturbation at the moment - adapt 
                                // for others later). "random" creates random initial conditions. "boost" creates a single straight (z directed) string with a Lorentz boost applied.
                                // "loop collision" creates two sets of (seperated) string, anti-string pairs. They are boosted towards each other and patched together so that they will collide
                                // and form a loop (2 one due to periodic boundary conditions which are required for this sim) 

const int nx = 201;
const int ny = 201;
const int nz = 51;
const double dx = 0.7;
const double dy = 0.7;
const double dz = 0.7;

const int n = 1; // This is useless for now, code assumes it is 1.
const double g = 1;

const double pi = 4*atan(1);

// Needed for straight string solutions based initial conditions /////////////////////

const int SORnx = 20001;
const double SORa = 0.01;

// Needed for perturbed straight strings (i.e sine string) ///////////////////////////

const double tol = 1e-6;

//const double eps = 0.1;
//const double lambda = dz*(nz-1)/(1-0.25*eps*eps);

const double gr = 0.5*(1+sqrt(5));

// Needed for random initial conditions ///////////////////////////////////////////////

const double seed = 1;

// Needed for boosted strings ///////////////////////////////////////////////////////////////////

const double dt = 0.2;

// Needed for the case of one z directed boosted string

// Relativistic velocities (c=1) so vx^2+vy^2+vz^2 must be less than 1

const double v1x = 0.6*0.4;
const double v1y = 0.6*sqrt(1-pow(0.4,2));
const double v1z = 0;

// Needed for loop collision. Assumes string anti-string pair 1 are x directed strings and pair 2 are z directed strings.

const double pos1s[2] = {10*0.6*sqrt(1-pow(0.4,2)), -25}; // y and z coordinates
const double pos1a[2] = {-10*0.6*sqrt(1-pow(0.4,2)), 25};

const double pos2s[2] = {25, -10*0.33*sqrt(1-pow(0.4,2))}; // x and y coordinates
const double pos2a[2] = {-25, 10*0.33*sqrt(1-pow(0.4,2))}; 


const double v1s[3] = {0, -0.6*sqrt(1-pow(0.4,2)), 0.6*0.4};
const double v1a[3] = {0, 0.6*sqrt(1-pow(0.4,2)), -0.6*0.4};

const double v2s[3] = {-0.33*0.4, 0.33*sqrt(1-pow(0.4,2)), 0};
const double v2a[3] = {0.33*0.4, -0.33*sqrt(1-pow(0.4,2)), 0};

const double omega = 0.5; // Phase modification parameter. Phase goes to zero more quickly for larger values
const double Lmod = 25; // Phase modification parameter. Length scale associated with modification


int main(){

    // User input to define epsilon

    double eps, lambda;
    cin >> eps;
    lambda = dz*(nz-1)/(1-0.25*eps*eps);

    Array phi(2,nx,ny,nz,0.0), A(3,nx,ny,nz,0.0);
    int i,j,k,comp;

    struct timeval start, end;
    gettimeofday(&start, NULL);

    string file_path = __FILE__;
    string dir_path = file_path.substr(0,file_path.find_last_of('/'));

    string icPath = dir_path + "/Data/ic.txt";
    string test1sPath = dir_path + "/test1s.txt";
    string test1aPath = dir_path + "/test1a.txt";
    string test2sPath = dir_path + "/test2s.txt";
    string test2aPath = dir_path + "/test2a.txt";

    ofstream ic (icPath.c_str());
    ofstream test1s (test1sPath.c_str());
    ofstream test1a (test1aPath.c_str());
    ofstream test2s (test2sPath.c_str());
    ofstream test2a (test2aPath.c_str());

    if(ic_type == "string"){

    	Array SOR_Fields(SORnx,2,0.0);

        double sigma[nx][nz], Fsigma[2], x, y, z, tolTest, Omega, distance, phiMag, AMag, normal_dist, xs, zs, xs_sigma, zs_sigma, xs_sigma2, zs_sigma2,
               paraVecMag, paraVecMag_sigma, a, b, sigma1, sigma2, distanceSqr1, distanceSqr2,finaldist;
        int pClosest;

        string SOR_inputPath = dir_path + "/SOR_Fields.txt";
        string testPath = dir_path + "/test.txt";

        ifstream SOR_input (SOR_inputPath.c_str());
        ofstream test (testPath.c_str());

        for(i=0;i<SORnx;i++){

    	   SOR_input >> SOR_Fields(i,0) >> SOR_Fields(i,1);

        }

        // Define epsilon and Omega based off amplitude and wavelength. Make sure that epsilon is small for accuracy. Maybe add check that it's at least less than 1.

        Omega = 2*pi/lambda;

        //cout << 2*eps/Omega << endl;

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

                // +1 in x distance so that z boundaries are fixed at x=0

                distanceSqr1 = pow(x - eps*(1 + cos(Omega*sigma1))/Omega,2) + pow(z - ( (1-0.25*eps*eps)*Omega*sigma1 + 0.125*eps*eps*sin(2*Omega*sigma1) )/Omega,2);
                distanceSqr2 = pow(x - eps*(1 + cos(Omega*sigma2))/Omega,2) + pow(z - ( (1-0.25*eps*eps)*Omega*sigma2 + 0.125*eps*eps*sin(2*Omega*sigma2) )/Omega,2);

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
   					
                        distanceSqr1 = pow(x - eps*(1 +cos(Omega*sigma1))/Omega,2) + pow(z - ( (1-0.25*eps*eps)*Omega*sigma1 + 0.125*eps*eps*sin(2*Omega*sigma1) )/Omega,2);
                        distanceSqr2 = pow(x - eps*(1 +cos(Omega*sigma2))/Omega,2) + pow(z - ( (1-0.25*eps*eps)*Omega*sigma2 + 0.125*eps*eps*sin(2*Omega*sigma2) )/Omega,2);

                    } else{

   					    a = sigma1;

   				        //cout << "2>1 and " << b - gr*(b-a) << " " << sigma1 << endl;

   				        sigma1 = b - (b-a)/gr;
   			      		sigma2 = a + (b-a)/gr;

   				     	distanceSqr1 = pow(x - eps*(1 +cos(Omega*sigma1))/Omega,2) + pow(z - ( (1-0.25*eps*eps)*Omega*sigma1 + 0.125*eps*eps*sin(2*Omega*sigma1) )/Omega,2);
   					    distanceSqr2 = pow(x - eps*(1 +cos(Omega*sigma2))/Omega,2) + pow(z - ( (1-0.25*eps*eps)*Omega*sigma2 + 0.125*eps*eps*sin(2*Omega*sigma2) )/Omega,2);

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

        		xs = eps*(1 + cos(Omega*sigma[i][k]))/Omega;
   		       	zs = ( (1-0.25*eps*eps)*Omega*sigma[i][k] + 0.125*eps*eps*sin(2*Omega*sigma[i][k]) )/Omega;

   	    		xs_sigma = -eps*sin(Omega*sigma[i][k]);
   	    		zs_sigma = 1-0.25*eps*eps + 0.25*eps*eps*cos(2*Omega*sigma[i][k]);

	       		paraVecMag = sqrt(pow(xs_sigma,2) + pow(zs_sigma,2));

    	       	normal_dist = ( zs_sigma*(x-xs) - xs_sigma*(z-zs) )/paraVecMag;

        		test << x << " " << z << " " << normal_dist << " " << sqrt(pow(x-xs,2) + pow(z-zs,2)) << " " << sqrt(finaldist) << endl;


    	       	for(j=0;j<ny;j++){

    	       		y = (j-y0)*dy;

    	       		distance = sqrt(pow(x-xs,2) + pow(y,2) + pow(z-zs,2));

    	       		pClosest = round(distance/SORa);

    		      	if(pClosest==0){

    		      		// 1st order interpolation since only have grid points on one side.

     	      			phiMag = ( SOR_Fields(pClosest+1,0)*distance - SOR_Fields(pClosest,0)*(distance-SORa) )/SORa;
     	      			AMag = ( SOR_Fields(pClosest+1,1)*distance - SOR_Fields(pClosest,1)*(distance-SORa) )/SORa;

     	      		} else if(pClosest<SORnx){

     	      			// 2nd order interpolation

     	      			phiMag = ( SOR_Fields(pClosest-1,0)*(distance-pClosest*SORa)*(distance-(pClosest+1)*SORa) - 2*SOR_Fields(pClosest,0)*(distance-(pClosest-1)*SORa)*(distance-(pClosest+1)*SORa) +
     		     				   SOR_Fields(pClosest+1,0)*(distance-(pClosest-1)*SORa)*(distance-pClosest*SORa) )/(2*SORa*SORa);

     		     		AMag = ( SOR_Fields(pClosest-1,1)*(distance-pClosest*SORa)*(distance-(pClosest+1)*SORa) - 2*SOR_Fields(pClosest,1)*(distance-(pClosest-1)*SORa)*(distance-(pClosest+1)*SORa) +
     		     				 SOR_Fields(pClosest+1,1)*(distance-(pClosest-1)*SORa)*(distance-pClosest*SORa) )/(2*SORa*SORa);

     	      		} else{

        				phiMag = 1;

     	      			if(g==0){ AMag = 0; }
        				else{ AMag = n/g; }

     	      			cout << "Off straight string solution grid" << endl;

        			}

        			// Now need to set phase of phi and split A_theta (theta defined with respect to string) into cartesian components


        			if(normal_dist==0){

     	      			// To prevent division by zero

     			    	phi(0,i,j,k) = 0;
     			    	A(1,i,j,k) = 0;

     			    } else{

     			    	phi(0,i,j,k) = phiMag*normal_dist/distance;
     		     		A(1,i,j,k) = AMag*normal_dist/pow(distance,2);

     	      		}

     	      		if(y==0){

     		     		// To prevent division by zero

     		     		phi(1,i,j,k) = 0;
     		     		A(0,i,j,k) = 0;
     	      			A(2,i,j,k) = 0;

     	      		} else{

     	      			phi(1,i,j,k) = phiMag*y/distance;
        				A(0,i,j,k) = -AMag*zs_sigma*y/(pow(distance,2)*paraVecMag);
     	      			A(2,i,j,k) = AMag*xs_sigma*y/(pow(distance,2)*paraVecMag);

     	      		}

     	      		// As long as this is not the first point in z (due to periodic boundary conditions) assign to z>0 using reflection symmetry across z=0

     	      		// if(k!=0){

     	      		// 	phi(0,i,j,nz-k) = phi(0,i,j,k);
        			// 	phi(1,i,j,nz-k) = phi(1,i,j,k);

        			// }

     	      		// Reflection across z=0

     	      		phi(0,i,j,nz-1-k) = phi(0,i,j,k);
     	      		phi(1,i,j,nz-1-k) = phi(1,i,j,k);

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

                ic << phi(0,i,j,k) << " " << phi(1,i,j,k) << " " << dx*g*A(0,i,j,k) << " " << dy*g*A(1,i,j,k) << " " << dz*g*A(2,i,j,k) << endl;

            }
        }
    }

    } else if(ic_type == "random"){

        uniform_real_distribution<double> unif_rand(-1,1);
        minstd_rand eng;

        eng.seed(seed); // Creates (psuedo) random number generator unif_rand with the seed specified at the beginning. Each number generated is between -1 and 1

        for(i=0;i<nx;i++){
            for(j=0;j<ny;j++){
                for(k=0;k<nz;k++){

                    phi(0,i,j,k) = unif_rand(eng);
                    phi(1,i,j,k) = unif_rand(eng);

                    if(g != 0){

                        // Gauge fields scaled by 1/g to try to account for "natural magnitude" of field. May need to adjust the way this is done

                        A(0,i,j,k) = unif_rand(eng)/g;
                        A(1,i,j,k) = unif_rand(eng)/g;
                        A(2,i,j,k) = unif_rand(eng)/g;

                    }

                    // Otherwise leave the A array as initialised - zeros everywhere

                    ic << phi(0,i,j,k) << " " << phi(1,i,j,k) << " " << dx*g*A(0,i,j,k) << " " << dy*g*A(1,i,j,k) << " " << dz*g*A(2,i,j,k) << endl;

                }
            }
        }

    } else if(ic_type == "boost"){

        Array SOR_Fields(SORnx,2,0.0), At(nx,ny,nz,0.0);

        double distance, phiMag, AMag, xb, yb, zb, xs, ys, zs, Ax, Ay, Re_phi, Im_phi;
        int pClosest;

        string SOR_inputPath = dir_path + "/SOR_Fields.txt";

        ifstream SOR_input (SOR_inputPath.c_str());

        for(i=0;i<SORnx;i++){

           SOR_input >> SOR_Fields(i,0) >> SOR_Fields(i,1);

        }

        // Grid positions of zero.

        int x0 = round((nx-1)/2);
        int y0 = round((ny-1)/2);
        int z0 = round((nz-1)/2);

        double vMagSqr = pow(v1x,2) + pow(v1y,2) + pow(v1z,2);

        if(vMagSqr>=1){ cout << "Error: boost velocity is greater than the speed of light" << endl; }

        double gamma = 1/sqrt(1 - vMagSqr);

        for(i=0;i<nx;i++){

            xb = (i-x0)*dx;

            for(j=0;j<ny;j++){

                yb = (j-y0)*dy;

                for(k=0;k<nz;k++){

                    zb = (k-z0)*dy;


                    // Firstly need to calculate what these positions correspond to in the frame with a stationary string at t=0. Don't care about z coordinate as symmetric

                    xs = xb + (gamma-1)*(v1x*xb + v1y*yb + v1z*zb)*v1x/vMagSqr;
                    ys = yb + (gamma-1)*(v1x*xb + v1y*yb + v1z*zb)*v1y/vMagSqr;

                    // Now calculate the distance from the string and interpolate the fields

                    distance = sqrt(pow(xs,2) + pow(ys,2));

                    pClosest = round(distance/SORa);

                    if(pClosest==0){

                        // 1st order interpolation since only have grid points on one side.

                        phiMag = ( SOR_Fields(pClosest+1,0)*distance - SOR_Fields(pClosest,0)*(distance-SORa) )/SORa;
                        AMag = ( SOR_Fields(pClosest+1,1)*distance - SOR_Fields(pClosest,1)*(distance-SORa) )/SORa;

                    } else if(pClosest<SORnx){

                        // 2nd order interpolation

                        phiMag = ( SOR_Fields(pClosest-1,0)*(distance-pClosest*SORa)*(distance-(pClosest+1)*SORa) - 2*SOR_Fields(pClosest,0)*(distance-(pClosest-1)*SORa)*(distance-(pClosest+1)*SORa) +
                                   SOR_Fields(pClosest+1,0)*(distance-(pClosest-1)*SORa)*(distance-pClosest*SORa) )/(2*SORa*SORa);

                        AMag = ( SOR_Fields(pClosest-1,1)*(distance-pClosest*SORa)*(distance-(pClosest+1)*SORa) - 2*SOR_Fields(pClosest,1)*(distance-(pClosest-1)*SORa)*(distance-(pClosest+1)*SORa) +
                                 SOR_Fields(pClosest+1,1)*(distance-(pClosest-1)*SORa)*(distance-pClosest*SORa) )/(2*SORa*SORa);

                    } else{

                        phiMag = 1;

                        if(g==0){ AMag = 0; }
                        else{ AMag = n/g; }

                        cout << "Off straight string solution grid" << endl;

                    }

                    // Now need to set phase of phi and split A_theta (theta defined with respect to string) into cartesian components. Az is zero everywhere for a z directed string so removed.


                    if(xs==0){

                        // To prevent division by zero

                        phi(0,i,j,k) = 0;
                        Ay = 0;

                    } else{

                        phi(0,i,j,k) = phiMag*xs/distance;
                        Ay = AMag*xs/pow(distance,2);

                    }

                    if(ys==0){

                        // To prevent division by zero

                        phi(1,i,j,k) = 0;
                        Ax = 0;

                    } else{

                        phi(1,i,j,k) = phiMag*ys/distance;
                        Ax = -AMag*ys/pow(distance,2);

                    }

                    // Now need to transform the gauge fields as one-forms under the Lorentz transformation. It is the inverse of the coordinate transformation but no difference at t=0.

                    At(i,j,k) = -gamma*(v1x*Ax + v1y*Ay);
                    A(0,i,j,k) = Ax + (gamma-1)*(v1x*Ax + v1y*Ay)*v1x/vMagSqr;
                    A(1,i,j,k) = Ay + (gamma-1)*(v1x*Ax + v1y*Ay)*v1y/vMagSqr;
                    A(2,i,j,k) = (gamma-1)*(v1x*Ax + v1y*Ay)*v1z/vMagSqr;

                    // At is needed for the gauge transformation at the next time step. Gauge transformation sets At(t=0) to zero without affecting Ai(t=0) but will have an effect on Ai(t=dt)

                    ic << phi(0,i,j,k) << " " << phi(1,i,j,k) << " " << dx*g*A(0,i,j,k) << " " << dy*g*A(1,i,j,k) << " " << dz*g*A(2,i,j,k) << endl;

                }
            }
        }


        // Now need to calculate the fields at the next timestep. This is slightly more complicated to do and also requires a gauge transformation to remain in the temporal gauge.
        // This is the reason for a separate loop as need to have calculated all At first.

        for(i=0;i<nx;i++){

            xb = (i-x0)*dx;

            for(j=0;j<ny;j++){

                yb = (j-y0)*dy;

                for(k=0;k<nz;k++){

                    zb = (k-z0)*dz;


                    // First calculate the positions in the stationary frame by Lorentz transforming from the moving frame. 
                    // xb - x coordinates in boosted frame. xs - x coordinates in stationary frame

                    xs = xb + (gamma-1)*(v1x*xb + v1y*yb + v1z*zb)*v1x/vMagSqr - gamma*v1x*dt;
                    ys = yb + (gamma-1)*(v1x*xb + v1y*yb + v1z*zb)*v1y/vMagSqr - gamma*v1y*dt;

                    // Now calculate the distance from the string and interpolate the fields

                    distance = sqrt(pow(xs,2) + pow(ys,2));

                    pClosest = round(distance/SORa);

                    if(pClosest==0){

                        // 1st order interpolation since only have grid points on one side.

                        phiMag = ( SOR_Fields(pClosest+1,0)*distance - SOR_Fields(pClosest,0)*(distance-SORa) )/SORa;
                        AMag = ( SOR_Fields(pClosest+1,1)*distance - SOR_Fields(pClosest,1)*(distance-SORa) )/SORa;

                    } else if(pClosest<SORnx){

                        // 2nd order interpolation

                        phiMag = ( SOR_Fields(pClosest-1,0)*(distance-pClosest*SORa)*(distance-(pClosest+1)*SORa) - 2*SOR_Fields(pClosest,0)*(distance-(pClosest-1)*SORa)*(distance-(pClosest+1)*SORa) +
                                   SOR_Fields(pClosest+1,0)*(distance-(pClosest-1)*SORa)*(distance-pClosest*SORa) )/(2*SORa*SORa);

                        AMag = ( SOR_Fields(pClosest-1,1)*(distance-pClosest*SORa)*(distance-(pClosest+1)*SORa) - 2*SOR_Fields(pClosest,1)*(distance-(pClosest-1)*SORa)*(distance-(pClosest+1)*SORa) +
                                 SOR_Fields(pClosest+1,1)*(distance-(pClosest-1)*SORa)*(distance-pClosest*SORa) )/(2*SORa*SORa);

                    } else{

                        phiMag = 1;

                        if(g==0){ AMag = 0; }
                        else{ AMag = n/g; }

                        cout << "Off straight string solution grid after initial timestep" << endl;

                    }

                    // Now need to set phase of phi and split A_theta (theta defined with respect to string) into cartesian components. Az is zero everywhere for a z directed string so removed.


                    if(xs==0){

                        // To prevent division by zero

                        Re_phi = 0;
                        Ay = 0;

                    } else{

                        Re_phi = phiMag*xs/distance;
                        Ay = AMag*xs/pow(distance,2);

                    }

                    if(ys==0){

                        // To prevent division by zero

                        Im_phi = 0;
                        Ax = 0;

                    } else{

                        Im_phi = phiMag*ys/distance;
                        Ax = -AMag*ys/pow(distance,2);

                    }

                    // Now need to transform the gauge fields as one-forms under the Lorentz transformation. Inverse of static to moving so same as coordinate transformations above.
                    // Since there is no time component of the gauge field in stationary frame the time part does not contribute unlike for the coordinates


                    A(0,i,j,k) = Ax + (gamma-1)*(v1x*Ax + v1y*Ay)*v1x/vMagSqr;
                    A(1,i,j,k) = Ay + (gamma-1)*(v1x*Ax + v1y*Ay)*v1y/vMagSqr;
                    A(2,i,j,k) = (gamma-1)*(v1x*Ax + v1y*Ay)*v1z/vMagSqr;

                    // Need to make a gauge transformation to remain in the temporal gauge after the gauge transformation. Assumes there are not periodic boundary conditions (as just one string)

                    // x component

                    if(i==0){

                        A(0,i,j,k) += -dt*( At(i+1,j,k) - At(i,j,k) )/dx;

                    } else if(i==nx-1){

                        A(0,i,j,k) += -dt*( At(i,j,k) - At(i-1,j,k) )/dx;

                    } else{

                        A(0,i,j,k) += -dt*( At(i+1,j,k) - At(i-1,j,k) )/(2*dx);

                    }

                    // y component

                    if(j==0){

                        A(1,i,j,k) += -dt*( At(i,j+1,k) - At(i,j,k) )/dy;

                    } else if(j==ny-1){

                        A(1,i,j,k) += -dt*( At(i,j,k) - At(i,j-1,k) )/dy;

                    } else{

                        A(1,i,j,k) += -dt*( At(i,j+1,k) - At(i,j-1,k) )/(2*dy);

                    }

                    // z component

                    if(k==0){

                        A(2,i,j,k) += -dt*( At(i,j,k+1) - At(i,j,k) )/dz;

                    } else if(k==nz-1){

                        A(2,i,j,k) += -dt*( At(i,j,k) - At(i,j,k-1) )/dz;

                    } else{

                        A(2,i,j,k) += -dt*( At(i,j,k+1) - At(i,j,k-1) )/(2*dz);

                    }

                    // Now transform the scalar field as well

                    phi(0,i,j,k) = cos(-g*dt*At(i,j,k))*Re_phi - sin(-g*dt*At(i,j,k))*Im_phi;
                    phi(1,i,j,k) = cos(-g*dt*At(i,j,k))*Im_phi + sin(-g*dt*At(i,j,k))*Re_phi;

                    ic << phi(0,i,j,k) << " " << phi(1,i,j,k) << " " << dx*g*A(0,i,j,k) << " " << dy*g*A(1,i,j,k) << " " << dz*g*A(2,i,j,k) << endl;

                }
            }
        }    

    } else if(ic_type == "loop collision"){

        // Could probably simplify a lot of this by incorporating all strings into one array. i.e make At(4,nx,ny,nz,0.0)

        Array SOR_Fields(SORnx,2,0.0), At1s(nx,ny,nz,0.0), At1a(nx,ny,nz,0.0), At2s(nx,ny,nz,0.0), At2a(nx,ny,nz,0.0);

        double distance[4], phiMag[4], AMag[4], xb, yb, zb, xs2s, xs2a, ys1s, ys1a, ys2s, ys2a, zs1s, zs1a, A1s[2], A1a[3], A2s[3], A2a[3], phi1s[2], phi1a[2], phi2s[2], phi2a[2],
               A1sPatch[3], A1aPatch[3], A2sPatch[3], A2aPatch[3];
        int pClosest[4], str;
        dcmplx ci(0.0,1.0), cUnitPos1s, cUnitPos1a, cUnitPos2s, cUnitPos2a, gRot1s, gRot1a, gRot2s, gRot2a;

        string SOR_inputPath = dir_path + "/SOR_Fields.txt";

        ifstream SOR_input (SOR_inputPath.c_str());

        for(i=0;i<SORnx;i++){

           SOR_input >> SOR_Fields(i,0) >> SOR_Fields(i,1);

        }

        // Grid positions of zero.

        int x0 = round((nx-1)/2);
        int y0 = round((ny-1)/2);
        int z0 = round((nz-1)/2);

        double v1sMagSqr = pow(v1s[0],2) + pow(v1s[1],2) + pow(v1s[2],2);
        double v1aMagSqr = pow(v1a[0],2) + pow(v1a[1],2) + pow(v1a[2],2);
        double v2sMagSqr = pow(v2s[0],2) + pow(v2s[1],2) + pow(v2s[2],2);
        double v2aMagSqr = pow(v2a[0],2) + pow(v2a[1],2) + pow(v2a[2],2);

        if(v1sMagSqr>=1 or v1aMagSqr>=1 or v2sMagSqr>=1 or v2aMagSqr>=1){ cout << "Error: boost velocity is greater than the speed of light" << endl; }

        double gamma1s = 1/sqrt(1 - v1sMagSqr);
        double gamma1a = 1/sqrt(1 - v1aMagSqr);
        double gamma2s = 1/sqrt(1 - v2sMagSqr);
        double gamma2a = 1/sqrt(1 - v2aMagSqr);

        for(i=0;i<nx;i++){

            xb = (i-x0)*dx;

            for(j=0;j<ny;j++){

                yb = (j-y0)*dy;

                for(k=0;k<nz;k++){

                    zb = (k-z0)*dy;


                    // Firstly need to calculate what these positions correspond to in the frame with a stationary string at t=0. Don't care about z coordinate as symmetric

                    xs2s = xb + (gamma2s-1)*(v2s[0]*xb + v2s[1]*yb + v2s[2]*zb)*v2s[0]/v2sMagSqr;
                    xs2a = xb + (gamma2a-1)*(v2a[0]*xb + v2a[1]*yb + v2a[2]*zb)*v2a[1]/v2aMagSqr;

                    ys1s = yb + (gamma1s-1)*(v1s[0]*xb - v1s[1]*yb + v1s[2]*zb)*v1s[1]/v1sMagSqr;
                    ys1a = yb + (gamma1a-1)*(v1a[0]*xb + v1a[1]*yb + v1a[2]*zb)*v1a[1]/v1aMagSqr;
                    ys2s = yb + (gamma2s-1)*(v2s[0]*xb + v2s[1]*yb + v2s[2]*zb)*v2s[1]/v2sMagSqr;
                    ys2a = yb - (gamma2a-1)*(v2a[0]*xb - v2a[1]*yb - v2a[2]*zb)*v2a[1]/v2aMagSqr;

                    zs1s = zb + (gamma1s-1)*(v1s[0]*xb - v1s[1]*yb + v1s[2]*zb)*v1s[2]/v1sMagSqr;
                    zs1a = zb + (gamma1a-1)*(v1a[0]*xb - v1a[1]*yb + v1a[2]*zb)*v1a[2]/v1aMagSqr;

                    // Now calculate the distance from the string and interpolate the fields

                    distance[0] = sqrt(pow(ys1s - pos1s[0],2) + pow(zs1s - pos1s[1],2));
                    distance[1] = sqrt(pow(ys1a - pos1a[0],2) + pow(zs1a - pos1a[1],2)); 
                    distance[2] = sqrt(pow(xs2s - pos2s[0],2) + pow(ys2s - pos2s[1],2)); 
                    distance[3] = sqrt(pow(xs2a - pos2a[0],2) + pow(ys2a - pos2a[1],2)); 

                    for(str=0;str<4;str++){

                        pClosest[str] = round(distance[str]/SORa);

                        if(pClosest[str]==0){

                            // 1st order interpolation since only have grid points on one side.

                            phiMag[str] = ( SOR_Fields(pClosest[str]+1,0)*distance[str] - SOR_Fields(pClosest[str],0)*(distance[str]-SORa) )/SORa;
                            AMag[str] = ( SOR_Fields(pClosest[str]+1,1)*distance[str] - SOR_Fields(pClosest[str],1)*(distance[str]-SORa) )/SORa;

                        } else if(pClosest[str]<SORnx){

                            // 2nd order interpolation

                            phiMag[str] = ( SOR_Fields(pClosest[str]-1,0)*(distance[str]-pClosest[str]*SORa)*(distance[str]-(pClosest[str]+1)*SORa) -
                                            2*SOR_Fields(pClosest[str],0)*(distance[str]-(pClosest[str]-1)*SORa)*(distance[str]-(pClosest[str]+1)*SORa) +
                                            SOR_Fields(pClosest[str]+1,0)*(distance[str]-(pClosest[str]-1)*SORa)*(distance[str]-pClosest[str]*SORa) )/(2*SORa*SORa);

                            AMag[str] = ( SOR_Fields(pClosest[str]-1,1)*(distance[str]-pClosest[str]*SORa)*(distance[str]-(pClosest[str]+1)*SORa) - 
                                          2*SOR_Fields(pClosest[str],1)*(distance[str]-(pClosest[str]-1)*SORa)*(distance[str]-(pClosest[str]+1)*SORa) +
                                          SOR_Fields(pClosest[str]+1,1)*(distance[str]-(pClosest[str]-1)*SORa)*(distance[str]-pClosest[str]*SORa) )/(2*SORa*SORa);

                        } else{

                            phiMag[str] = 1;

                            if(g==0){ AMag[str] = 0; }
                            else{ AMag[str] = n/g; }

                            cout << "Off straight string solution grid" << endl;

                        }

                    }

                    // Now need to set phase of phi and split A_theta (theta defined with respect to string) into cartesian components for each string. Modify phase for patching at the same time

                    cUnitPos1s = ( -(zs1s - pos1s[1]) + ci*(ys1s - pos1s[0]) )/distance[0];
                    cUnitPos1a = ( (zs1a - pos1a[1]) + ci*(ys1a - pos1a[0]) )/distance[1];
                    cUnitPos2s = ( (xs2s - pos2s[0]) + ci*(ys2s - pos2s[1]) )/distance[2];
                    cUnitPos2a = ( -(xs2a - pos2a[0]) + ci*(ys2a - pos2a[1]) )/distance[3];

                    // string 1

                    if(zs1s - pos1s[1] == 0){ phi1s[0] = 0; A1s[1] = 0; } // To prevent division by zero
                    else{

                        phi1s[0] = phiMag[0]*real( pow(cUnitPos1s, 0.5*( 1 - tanh(omega*(distance[0] - Lmod)) )) ); // The pow() is the modification for patching
                        A1s[1] = -AMag[0]*(zs1s - pos1s[1])/pow(distance[0],2); // y component

                    }

                    if(ys1s - pos1s[0] == 0){ phi1s[1] = 0; A1s[2] = 0; }
                    else{

                        phi1s[1] = phiMag[0]*imag( pow(cUnitPos1s, 0.5*( 1 - tanh(omega*(distance[0] - Lmod)) )) );
                        A1s[2] = AMag[0]*(ys1s - pos1s[0])/pow(distance[0],2); // z component

                    }

                    // antistring 1

                    if(zs1a - pos1a[1] == 0){ phi1a[0] = 0; A1a[1] = 0; }
                    else{

                        phi1a[0] = phiMag[1]*real( pow(cUnitPos1a, 0.5*( 1 - tanh(omega*(distance[1] - Lmod)) )) );
                        A1a[1] = AMag[1]*(zs1a - pos1a[1])/pow(distance[1],2); // y component

                    }

                    if(ys1a - pos1a[0] == 0){ phi1a[1] = 0; A1a[2] = 0; }
                    else{

                        phi1a[1] = phiMag[1]*imag( pow(cUnitPos1a, 0.5*( 1 - tanh(omega*(distance[1] - Lmod)) )) );
                        A1a[2] = -AMag[1]*(ys1a - pos1a[0])/pow(distance[1],2); // z component

                    }

                    // string 2

                    if(xs2s - pos2s[0] == 0){ phi2s[0] = 0; A2s[1] = 0; }
                    else{

                        phi2s[0] = phiMag[2]*real( pow(cUnitPos2s, 0.5*( 1 - tanh(omega*(distance[2] - Lmod)) )) );
                        A2s[1] = AMag[2]*(xs2s - pos2s[0])/pow(distance[2],2); // y component

                    }

                    if(ys2s - pos2s[1] == 0){ phi2s[1] = 0; A2s[0] = 0; }
                    else{

                        phi2s[1] = phiMag[2]*imag( pow(cUnitPos2s, 0.5*( 1 - tanh(omega*(distance[2] - Lmod)) )) );
                        A2s[0] = -AMag[2]*(ys2s - pos2s[1])/pow(distance[2],2); // x component

                    }

                    // antistring 2

                    if(xs2a - pos2a[0] == 0){ phi2a[0] = 0; A2a[1] = 0; }
                    else{

                        phi2a[0] = phiMag[3]*real( pow(cUnitPos2a, 0.5*( 1 - tanh(omega*(distance[3] - Lmod)) )) );
                        A2a[1] = -AMag[3]*(xs2a - pos2a[0])/pow(distance[3],2); // y component

                    }

                    if(ys2a - pos2a[1] == 0){ phi2a[1] = 0; A2a[0] = 0; }
                    else{

                        phi2a[1] = phiMag[3]*imag( pow(cUnitPos2a, 0.5*( 1 - tanh(omega*(distance[3] - Lmod)) )) );
                        A2a[0] = AMag[3]*(ys2a - pos2a[1])/pow(distance[3],2); // x component

                    }

                    // Now need to transform the gauge fields as one-forms under the Lorentz transformation. It is the inverse of the coordinate transformation but no difference at t=0.

                    At1s(i,j,k) = -gamma1s*(v1s[1]*A1s[1] + v1s[2]*A1s[2]);
                    At1a(i,j,k) = -gamma1a*(v1a[1]*A1a[1] + v1a[2]*A1a[2]);
                    At2s(i,j,k) = -gamma2s*(v2s[0]*A2s[0] + v2s[1]*A2s[1]);
                    At2a(i,j,k) = -gamma2a*(v2a[0]*A2a[0] + v2a[1]*A2a[1]);

                    // At is needed for the gauge transformation at the next time step. Gauge transformation sets At(t=0) to zero without affecting Ai(t=0) but will have an effect on Ai(t=dt)

                    for(comp=0;comp<3;comp++){ 

                        A1sPatch[comp] = (gamma1s-1)*(v1s[1]*A1s[1] + v1s[2]*A1s[2])*v1s[comp]/v1sMagSqr;
                        A1aPatch[comp] = (gamma1a-1)*(v1a[1]*A1a[1] + v1a[2]*A1a[2])*v1a[comp]/v1aMagSqr;
                        A2sPatch[comp] = (gamma2s-1)*(v2s[0]*A2s[0] + v2s[1]*A2s[1])*v2s[comp]/v2sMagSqr;
                        A2aPatch[comp] = (gamma2a-1)*(v2a[0]*A2a[0] + v2a[1]*A2a[1])*v2a[comp]/v2aMagSqr; 

                        // Patch the strings gauge fields together

                        A(comp,i,j,k) = A1sPatch[comp] + A1aPatch[comp] + A2sPatch[comp] + A2aPatch[comp];

                    }

                    // Patch the strings scalar fields together

                    phi(0,i,j,k) = phi1s[0]*phi1a[0]*phi2s[0]*phi2a[0] - phi1s[0]*phi1a[0]*phi2s[1]*phi2a[1] - phi1s[0]*phi1a[1]*phi2s[0]*phi2a[1] - phi1s[0]*phi1a[1]*phi2s[1]*phi2a[0]
                                 - phi1s[1]*phi1a[0]*phi2s[0]*phi2a[1] - phi1s[1]*phi1a[0]*phi2s[1]*phi2a[0] - phi1s[1]*phi1a[1]*phi2s[0]*phi2a[0] + phi1s[1]*phi1a[1]*phi2s[1]*phi2a[1];

                    phi(1,i,j,k) = phi1s[0]*phi1a[0]*phi2s[0]*phi2a[1] + phi1s[0]*phi1a[0]*phi2s[1]*phi2a[0] + phi1s[0]*phi1a[1]*phi2s[0]*phi2a[0] + phi1s[1]*phi1a[0]*phi2s[0]*phi2a[0]
                                 - phi1s[1]*phi1a[1]*phi2s[1]*phi2a[0] - phi1s[1]*phi1a[1]*phi2s[0]*phi2a[1] - phi1s[1]*phi1a[0]*phi2s[1]*phi2a[1] - phi1s[0]*phi1a[1]*phi2s[1]*phi2a[1];

                    ic << phi(0,i,j,k) << " " << phi(1,i,j,k) << " " << dx*g*A(0,i,j,k) << " " << dy*g*A(1,i,j,k) << " " << dz*g*A(2,i,j,k) << endl;

                }
            }
        }


        // Now set generate the field values at the next timestep


        for(i=0;i<nx;i++){

            xb = (i-x0)*dx;

            for(j=0;j<ny;j++){

                yb = (j-y0)*dy;

                for(k=0;k<nz;k++){

                    zb = (k-z0)*dy;


                    // Firstly need to calculate what these positions correspond to in the frame with a stationary string at t=0. Don't care about z coordinate as symmetric


                    xs2s = xb + (gamma2s-1)*(v2s[0]*xb + v2s[1]*yb + v2s[2]*zb)*v2s[0]/v2sMagSqr - gamma2s*v2s[0]*dt;
                    xs2a = xb + (gamma2a-1)*(v2a[0]*xb + v2a[1]*yb + v2a[2]*zb)*v2a[1]/v2aMagSqr - gamma2a*v2a[0]*dt;

                    ys1s = yb + (gamma1s-1)*(v1s[0]*xb - v1s[1]*yb + v1s[2]*zb)*v1s[1]/v1sMagSqr - gamma1s*v1s[1]*dt;
                    ys1a = yb + (gamma1a-1)*(v1a[0]*xb + v1a[1]*yb + v1a[2]*zb)*v1a[1]/v1aMagSqr - gamma1a*v1a[1]*dt;
                    ys2s = yb + (gamma2s-1)*(v2s[0]*xb + v2s[1]*yb + v2s[2]*zb)*v2s[1]/v2sMagSqr - gamma2s*v2s[1]*dt;
                    ys2a = yb - (gamma2a-1)*(v2a[0]*xb - v2a[1]*yb - v2a[2]*zb)*v2a[1]/v2aMagSqr - gamma2a*v2a[1]*dt;

                    zs1s = zb + (gamma1s-1)*(v1s[0]*xb - v1s[1]*yb + v1s[2]*zb)*v1s[2]/v1sMagSqr - gamma1s*v1s[2]*dt;
                    zs1a = zb + (gamma1a-1)*(v1a[0]*xb - v1a[1]*yb + v1a[2]*zb)*v1a[2]/v1aMagSqr - gamma1a*v1a[2]*dt;

                    // Now calculate the distance from the string and interpolate the fields

                    distance[0] = sqrt(pow(ys1s - pos1s[0],2) + pow(zs1s - pos1s[1],2));
                    distance[1] = sqrt(pow(ys1a - pos1a[0],2) + pow(zs1a - pos1a[1],2)); 
                    distance[2] = sqrt(pow(xs2s - pos2s[0],2) + pow(ys2s - pos2s[1],2)); 
                    distance[3] = sqrt(pow(xs2a - pos2a[0],2) + pow(ys2a - pos2a[1],2)); 

                    for(str=0;str<4;str++){

                        pClosest[str] = round(distance[str]/SORa);

                        if(pClosest[str]==0){

                            // 1st order interpolation since only have grid points on one side.

                            phiMag[str] = ( SOR_Fields(pClosest[str]+1,0)*distance[str] - SOR_Fields(pClosest[str],0)*(distance[str]-SORa) )/SORa;
                            AMag[str] = ( SOR_Fields(pClosest[str]+1,1)*distance[str] - SOR_Fields(pClosest[str],1)*(distance[str]-SORa) )/SORa;

                        } else if(pClosest[str]<SORnx){

                            // 2nd order interpolation

                            phiMag[str] = ( SOR_Fields(pClosest[str]-1,0)*(distance[str]-pClosest[str]*SORa)*(distance[str]-(pClosest[str]+1)*SORa) -
                                            2*SOR_Fields(pClosest[str],0)*(distance[str]-(pClosest[str]-1)*SORa)*(distance[str]-(pClosest[str]+1)*SORa) +
                                            SOR_Fields(pClosest[str]+1,0)*(distance[str]-(pClosest[str]-1)*SORa)*(distance[str]-pClosest[str]*SORa) )/(2*SORa*SORa);

                            AMag[str] = ( SOR_Fields(pClosest[str]-1,1)*(distance[str]-pClosest[str]*SORa)*(distance[str]-(pClosest[str]+1)*SORa) - 
                                          2*SOR_Fields(pClosest[str],1)*(distance[str]-(pClosest[str]-1)*SORa)*(distance[str]-(pClosest[str]+1)*SORa) +
                                          SOR_Fields(pClosest[str]+1,1)*(distance[str]-(pClosest[str]-1)*SORa)*(distance[str]-pClosest[str]*SORa) )/(2*SORa*SORa);

                        } else{

                            phiMag[str] = 1;

                            if(g==0){ AMag[str] = 0; }
                            else{ AMag[str] = n/g; }

                            cout << "Off straight string solution grid" << endl;

                        }

                    }

                    // Now need to set phase of phi and split A_theta (theta defined with respect to string) into cartesian components for each string. Modify phase for patching at the same time

                    cUnitPos1s = ( -(zs1s - pos1s[1]) + ci*(ys1s - pos1s[0]) )/distance[0];
                    cUnitPos1a = ( (zs1a - pos1a[1]) + ci*(ys1a - pos1a[0]) )/distance[1];
                    cUnitPos2s = ( (xs2s - pos2s[0]) + ci*(ys2s - pos2s[1]) )/distance[2];
                    cUnitPos2a = ( -(xs2a - pos2a[0]) + ci*(ys2a - pos2a[1]) )/distance[3];

                    // Gauge transformations of phi. phi --> phi*e^-i*g*dt*At(t=0)

                    gRot1s = cos(g*dt*At1s(i,j,k)) - ci*sin(g*dt*At1s(i,j,k));
                    gRot1a = cos(g*dt*At1a(i,j,k)) - ci*sin(g*dt*At1a(i,j,k));
                    gRot2s = cos(g*dt*At2s(i,j,k)) - ci*sin(g*dt*At2s(i,j,k));
                    gRot2a = cos(g*dt*At2a(i,j,k)) - ci*sin(g*dt*At2a(i,j,k));

                    // string 1

                    if(zs1s - pos1s[1] == 0){ phi1s[0] = 0; A1s[1] = 0; } // To prevent division by zero
                    else{

                        phi1s[0] = phiMag[0]*real( pow(cUnitPos1s*gRot1s, 0.5*( 1 - tanh(omega*(distance[0] - Lmod)) )) ); // The pow() is the modification for patching
                        A1s[1] = -AMag[0]*(zs1s - pos1s[1])/pow(distance[0],2); // y component

                    }

                    if(ys1s - pos1s[0] == 0){ phi1s[1] = 0; A1s[2] = 0; }
                    else{

                        phi1s[1] = phiMag[0]*imag( pow(cUnitPos1s*gRot1s, 0.5*( 1 - tanh(omega*(distance[0] - Lmod)) )) );
                        A1s[2] = AMag[0]*(ys1s - pos1s[0])/pow(distance[0],2); // z component

                    }

                    // antistring 1

                    if(zs1a - pos1a[1] == 0){ phi1a[0] = 0; A1a[1] = 0; }
                    else{

                        phi1a[0] = phiMag[1]*real( pow(cUnitPos1a*gRot1a, 0.5*( 1 - tanh(omega*(distance[1] - Lmod)) )) );
                        A1a[1] = AMag[1]*(zs1a - pos1a[1])/pow(distance[1],2); // y component

                    }

                    if(ys1a - pos1a[0] == 0){ phi1a[1] = 0; A1a[2] = 0; }
                    else{

                        phi1a[1] = phiMag[1]*imag( pow(cUnitPos1a*gRot1a, 0.5*( 1 - tanh(omega*(distance[1] - Lmod)) )) );
                        A1a[2] = -AMag[1]*(ys1a - pos1a[0])/pow(distance[1],2); // z component

                    }

                    // string 2

                    if(xs2s - pos2s[0] == 0){ phi2s[0] = 0; A2s[1] = 0; }
                    else{

                        phi2s[0] = phiMag[2]*real( pow(cUnitPos2s*gRot2s, 0.5*( 1 - tanh(omega*(distance[2] - Lmod)) )) );
                        A2s[1] = AMag[2]*(xs2s - pos2s[0])/pow(distance[2],2); // y component

                    }

                    if(ys2s - pos2s[1] == 0){ phi2s[1] = 0; A2s[0] = 0; }
                    else{

                        phi2s[1] = phiMag[2]*imag( pow(cUnitPos2s*gRot2s, 0.5*( 1 - tanh(omega*(distance[2] - Lmod)) )) );
                        A2s[0] = -AMag[2]*(ys2s - pos2s[1])/pow(distance[2],2); // x component

                    }

                    // antistring 2

                    if(xs2a - pos2a[0] == 0){ phi2a[0] = 0; A2a[1] = 0; }
                    else{

                        phi2a[0] = phiMag[3]*real( pow(cUnitPos2a*gRot2a, 0.5*( 1 - tanh(omega*(distance[3] - Lmod)) )) );
                        A2a[1] = -AMag[3]*(xs2a - pos2a[0])/pow(distance[3],2); // y component

                    }

                    if(ys2a - pos2a[1] == 0){ phi2a[1] = 0; A2a[0] = 0; }
                    else{

                        phi2a[1] = phiMag[3]*imag( pow(cUnitPos2a*gRot2a, 0.5*( 1 - tanh(omega*(distance[3] - Lmod)) )) );
                        A2a[0] = AMag[3]*(ys2a - pos2a[1])/pow(distance[3],2); // x component

                    }

                    // Lorentz transformation of the gauge field components

                    for(comp=0;comp<3;comp++){ 

                        A1sPatch[comp] = (gamma1s-1)*(v1s[1]*A1s[1] + v1s[2]*A1s[2])*v1s[comp]/v1sMagSqr;
                        A1aPatch[comp] = (gamma1a-1)*(v1a[1]*A1a[1] + v1a[2]*A1a[2])*v1a[comp]/v1aMagSqr;
                        A2sPatch[comp] = (gamma2s-1)*(v2s[0]*A2s[0] + v2s[1]*A2s[1])*v2s[comp]/v2sMagSqr;
                        A2aPatch[comp] = (gamma2a-1)*(v2a[0]*A2a[0] + v2a[1]*A2a[1])*v2a[comp]/v2aMagSqr; 

                    }

                    // Now need to make the gauge transformation of the gauge fields to set At=0. Assumes periodic boundary conditions

                    // x components

                    if(i==0){

                        A1sPatch[0] += -dt*( At1s(i+1,j,k) - At1s(nx-1,j,k) )/(2*dx);
                        A1aPatch[0] += -dt*( At1a(i+1,j,k) - At1a(nx-1,j,k) )/(2*dx);
                        A2sPatch[0] += -dt*( At2s(i+1,j,k) - At2s(nx-1,j,k) )/(2*dx);
                        A2aPatch[0] += -dt*( At2a(i+1,j,k) - At2a(nx-1,j,k) )/(2*dx);

                    } else if(i==nx-1){

                        A1sPatch[0] += -dt*( At1s(0,j,k) - At1s(i-1,j,k) )/(2*dx);
                        A1aPatch[0] += -dt*( At1a(0,j,k) - At1a(i-1,j,k) )/(2*dx);
                        A2sPatch[0] += -dt*( At2s(0,j,k) - At2s(i-1,j,k) )/(2*dx);
                        A2aPatch[0] += -dt*( At2a(0,j,k) - At2a(i-1,j,k) )/(2*dx);

                    } else{

                        A1sPatch[0] += -dt*( At1s(i+1,j,k) - At1s(i-1,j,k) )/(2*dx);
                        A1aPatch[0] += -dt*( At1a(i+1,j,k) - At1a(i-1,j,k) )/(2*dx);
                        A2sPatch[0] += -dt*( At2s(i+1,j,k) - At2s(i-1,j,k) )/(2*dx);
                        A2aPatch[0] += -dt*( At2a(i+1,j,k) - At2a(i-1,j,k) )/(2*dx);

                    }

                    // y components

                    if(j==0){

                        A1sPatch[1] += -dt*( At1s(i,j+1,k) - At1s(i,ny-1,k) )/(2*dy);
                        A1aPatch[1] += -dt*( At1a(i,j+1,k) - At1a(i,ny-1,k) )/(2*dy);
                        A2sPatch[1] += -dt*( At2s(i,j+1,k) - At2s(i,ny-1,k) )/(2*dy);
                        A2aPatch[1] += -dt*( At2a(i,j+1,k) - At2a(i,ny-1,k) )/(2*dy);

                    } else if(j==ny-1){

                        A1sPatch[1] += -dt*( At1s(i,0,k) - At1s(i,j-1,k) )/(2*dy);
                        A1aPatch[1] += -dt*( At1a(i,0,k) - At1a(i,j-1,k) )/(2*dy);
                        A2sPatch[1] += -dt*( At2s(i,0,k) - At2s(i,j-1,k) )/(2*dy);
                        A2aPatch[1] += -dt*( At2a(i,0,k) - At2a(i,j-1,k) )/(2*dy);

                    } else{

                        A1sPatch[1] += -dt*( At1s(i,j+1,k) - At1s(i,j-1,k) )/(2*dy);
                        A1aPatch[1] += -dt*( At1a(i,j+1,k) - At1a(i,j-1,k) )/(2*dy);
                        A2sPatch[1] += -dt*( At2s(i,j+1,k) - At2s(i,j-1,k) )/(2*dy);
                        A2aPatch[1] += -dt*( At2a(i,j+1,k) - At2a(i,j-1,k) )/(2*dy);

                    }

                    // z components

                    if(k==0){

                        A1sPatch[2] += -dt*( At1s(i,j,k+1) - At1s(i,j,nz-1) )/(2*dz);
                        A1aPatch[2] += -dt*( At1a(i,j,k+1) - At1a(i,j,nz-1) )/(2*dz);
                        A2sPatch[2] += -dt*( At2s(i,j,k+1) - At2s(i,j,nz-1) )/(2*dz);
                        A2aPatch[2] += -dt*( At2a(i,j,k+1) - At2a(i,j,nz-1) )/(2*dz);

                    } else if(k==nz-1){

                        A1sPatch[2] += -dt*( At1s(i,j,0) - At1s(i,j,k-1) )/(2*dz);
                        A1aPatch[2] += -dt*( At1a(i,j,0) - At1a(i,j,k-1) )/(2*dz);
                        A2sPatch[2] += -dt*( At2s(i,j,0) - At2s(i,j,k-1) )/(2*dz);
                        A2aPatch[2] += -dt*( At2a(i,j,0) - At2a(i,j,k-1) )/(2*dz);

                    } else{

                        A1sPatch[2] += -dt*( At1s(i,j,k+1) - At1s(i,j,k-1) )/(2*dz);
                        A1aPatch[2] += -dt*( At1a(i,j,k+1) - At1a(i,j,k-1) )/(2*dz);
                        A2sPatch[2] += -dt*( At2s(i,j,k+1) - At2s(i,j,k-1) )/(2*dz);
                        A2aPatch[2] += -dt*( At2a(i,j,k+1) - At2a(i,j,k-1) )/(2*dz);

                    }


                    // Patch the strings together

                    phi(0,i,j,k) = phi1s[0]*phi1a[0]*phi2s[0]*phi2a[0] - phi1s[0]*phi1a[0]*phi2s[1]*phi2a[1] - phi1s[0]*phi1a[1]*phi2s[0]*phi2a[1] - phi1s[0]*phi1a[1]*phi2s[1]*phi2a[0]
                                 - phi1s[1]*phi1a[0]*phi2s[0]*phi2a[1] - phi1s[1]*phi1a[0]*phi2s[1]*phi2a[0] - phi1s[1]*phi1a[1]*phi2s[0]*phi2a[0] + phi1s[1]*phi1a[1]*phi2s[1]*phi2a[1];

                    phi(1,i,j,k) = phi1s[0]*phi1a[0]*phi2s[0]*phi2a[1] + phi1s[0]*phi1a[0]*phi2s[1]*phi2a[0] + phi1s[0]*phi1a[1]*phi2s[0]*phi2a[0] + phi1s[1]*phi1a[0]*phi2s[0]*phi2a[0]
                                 - phi1s[1]*phi1a[1]*phi2s[1]*phi2a[0] - phi1s[1]*phi1a[1]*phi2s[0]*phi2a[1] - phi1s[1]*phi1a[0]*phi2s[1]*phi2a[1] - phi1s[0]*phi1a[1]*phi2s[1]*phi2a[1];

                    for(comp=0;comp<3;comp++){ A(comp,i,j,k) = A1sPatch[comp] + A1aPatch[comp] + A2sPatch[comp] + A2aPatch[comp]; }

                    ic << phi(0,i,j,k) << " " << phi(1,i,j,k) << " " << dx*g*A(0,i,j,k) << " " << dy*g*A(1,i,j,k) << " " << dz*g*A(2,i,j,k) << endl;

                    test1s << phi1s[0] << " " << phi1s[1] << endl;
                    test1a << phi1a[0] << " " << phi1a[1] << endl;
                    test2s << phi2s[0] << " " << phi2s[1] << endl;
                    test2a << phi2a[0] << " " << phi2a[1] << endl;

                }
            }
        }



    } else{

        cout << "Unrecognised ic_type" << endl;

    }


    gettimeofday(&end,NULL);

    cout << "Time taken: " << end.tv_sec - start.tv_sec << "s" << endl;


	return 0;

}


