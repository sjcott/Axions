#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <sys/time.h>
#include <omp.h>
#include <complex>
#include <vector>

#include "header.hpp"

using namespace std;

// Perturb a straight, global string solution, evolve the system and investigate how it radiates.

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                            Class Definitions
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

const int nx = 200;
const int ny = 200;
const int nz = 200;
const int nt = 6001;
const double dx = 0.25;
const double dy = 0.25;
const double dz = 0.25;
const double dt = 0.1;

const bool makeGif = false;
const int saveFreq = 60;


int main(){

	Array phi(2,nx,ny,nz,0.0), phit(2,nx,ny,nz,0.0), phitt(2,nx,ny,nz,0.0), energydensity(nx,ny,nz,0.0);
	int comp,i,j,k,TimeStep, gifFrame;
    double phixx, phiyy, phizz, phiMagSqr, phix[2], phiy[2], phiz[2], energy;

	struct timeval start, end;
    gettimeofday(&start, NULL);

    string file_path = __FILE__;
    string dir_path = file_path.substr(0,file_path.find_last_of('/'));
    stringstream ss;

    string icPath = dir_path + "/ic.txt";
    string finalFieldPath = dir_path + "/finalField.txt";
    string valsPerLoopPath = dir_path + "/valsPerLoop.txt";

    ifstream ic (icPath.c_str());
    ofstream finalField (finalFieldPath.c_str());
    ofstream valsPerLoop (valsPerLoopPath.c_str());

    for(i=0;i<nx;i++){
        for(j=0;j<ny;j++){
            for(k=0;k<nz;k++){

                ic >> phi(0,i,j,k) >> phi(1,i,j,k);

                // Time derivatives are initialised to zero by class constructor

            }
        }
    }

    // Calculate initial energy

    // energy = 0;

    // for(i=1;i<nx-1;i++){
    //     for(j=1;j<ny-1;j++){
    //         for(k=0;k<nz;k++){

    //             phiMagSqr = pow(phi(0,i,j,k),2) + pow(phi(1,i,j,k),2);

    //             for(comp=0;comp<2;comp++){

    //                 phix[comp] = ( phi(comp,i+1,j,k) - phi(comp,i-1,j,k) )/(2*a);
    //                 phiy[comp] = ( phi(comp,i,j+1,k) - phi(comp,i,j-1,k) )/(2*a);

    //                 if(k==0){

    //                     phiz[comp] = ( phi(comp,i,j,k+1) - phi(comp,i,j,nz-1) )/(2*a);

    //                 } else if(k==nz-1){

    //                     phiz[comp] = ( phi(comp,i,j,0) - phi(comp,i,j,k-1) )/(2*a);

    //                 } else{

    //                     phiz[comp] = ( phi(comp,i,j,k+1) - phi(comp,i,j,k-1) )/(2*a);

    //                 }

    //             }

    //             energydensity(i,j,k) = pow(phix[0],2) + pow(phix[1],2) + pow(phiy[0],2) + pow(phiy[1],2) + pow(phiz[0],2) + pow(phiz[1],2) + 0.5*pow(phiMagSqr-1,2);

    //             energy += pow(a,3)*energydensity(i,j,k);

    //         }
    //     }
    // }


    //valsPerLoop << energy << endl;

    gifFrame = 0;

    for(TimeStep=0;TimeStep<nt;TimeStep++){

        if((1000*TimeStep)%(nt-1)==0){

            cout << "\rTimestep " << TimeStep << " completed." << flush;

        }

        // Set boundary conditions. Different loop sizes so that corner allocation is done correctly.

        // Sets Neumann boundary conditions as simple first attempt. Will need to use absorbing boundary conditions in future. z direction has periodic boundary conditions instead.

        for(i=1;i<nx-1;i++){
            for(k=0;k<nz;k++){

                phi(0,i,0,k) = phi(0,i,1,k);
                phi(1,i,0,k) = phi(1,i,1,k);

                phi(0,i,ny-1,k) = phi(0,i,ny-2,k);
                phi(1,i,ny-1,k) = phi(1,i,ny-2,k);

            }
        }

        for(j=0;j<ny;j++){
            for(k=0;k<nz;k++){

                phi(0,0,j,k) = phi(0,1,j,k);
                phi(1,0,j,k) = phi(1,1,j,k);

                phi(0,nx-1,j,k) = phi(0,nx-2,j,k);
                phi(1,nx-1,j,k) = phi(1,nx-2,j,k);

            }
        }


        // Calculate time derivatives using EoMs

        #pragma omp parallel for default(none) shared(phitt,phi) private(phiMagSqr,phixx,phiyy,phizz,j,k,comp)
        for(i=1;i<nx-1;i++){
            for(j=1;j<ny-1;j++){
                for(k=0;k<nz;k++){

                    phiMagSqr = pow(phi(0,i,j,k),2) + pow(phi(1,i,j,k),2);

                    for(comp=0;comp<2;comp++){

                         // 2nd order spatial derivatives calculated with 2nd order finite difference

                        phixx = ( phi(comp,i+1,j,k) - 2*phi(comp,i,j,k) + phi(comp,i-1,j,k) )/(dx*dx);
                        phiyy = ( phi(comp,i,j+1,k) - 2*phi(comp,i,j,k) + phi(comp,i,j-1,k) )/(dy*dy);

                        // if(i==0){

                        //     phixx = ( phi(comp,i+1,j,k) - 2*phi(comp,i,j,k) + phi(comp,nx-1,j,k) )/(dx*dx);

                        // } else if(i==nx-1){

                        //     phixx = ( phi(comp,0,j,k) - 2*phi(comp,i,j,k) + phi(comp,i-1,j,k) )/(dx*dx);

                        // } else{

                        //     phixx = ( phi(comp,i+1,j,k) - 2*phi(comp,i,j,k) + phi(comp,i-1,j,k) )/(dx*dx);

                        // }

                        // if(j==0){

                        //     phiyy = ( phi(comp,i,j+1,k) - 2*phi(comp,i,j,k) + phi(comp,i,ny-1,k) )/(dy*dy);

                        // } else if(j==ny-1){

                        //     phiyy = ( phi(comp,i,0,k) - 2*phi(comp,i,j,k) + phi(comp,i,j-1,k) )/(dy*dy);

                        // } else{

                        //     phiyy = ( phi(comp,i,j+1,k) - 2*phi(comp,i,j,k) + phi(comp,i,j-1,k) )/(dy*dy);

                        // }

                        // enforce periodic boundary condition in z direction

                        if(k==0){

                            phizz = ( phi(comp,i,j,k+1) - 2*phi(comp,i,j,k) + phi(comp,i,j,nz-1) )/(dz*dz);

                        } else if(k==nz-1){

                            phizz = ( phi(comp,i,j,0) - 2*phi(comp,i,j,k) + phi(comp,i,j,k-1) )/(dz*dz);

                        } else{
                        
                            phizz = ( phi(comp,i,j,k+1) - 2*phi(comp,i,j,k) + phi(comp,i,j,k-1) )/(dz*dz);

                        }

                        // Calculate second order time derivatives

                        phitt(comp,i,j,k) = phixx + phiyy + phizz - (phiMagSqr - 1)*phi(comp,i,j,k);

                    }
                }
            }
        }


        // Update the arrays and calculate new energy

        energy = 0;

        #pragma omp parallel for reduction(+:energy) default(none) shared(energydensity,phi,phit,phitt) private(phiMagSqr,phix,phiy,phiz,j,k,comp)
        for(i=1;i<nx-1;i++){
            for(j=1;j<ny-1;j++){
                for(k=0;k<nz;k++){

                    phiMagSqr = pow(phi(0,i,j,k),2) + pow(phi(1,i,j,k),2);

                    for(comp=0;comp<2;comp++){

                        // Update momentum array

                        phit(comp,i,j,k) += phitt(comp,i,j,k)*dt;


                        // Calculate 1st derivatives for energy calculation

                        phix[comp] = ( phi(comp,i+1,j,k) - phi(comp,i-1,j,k) )/(2*dx);
                        phiy[comp] = ( phi(comp,i,j+1,k) - phi(comp,i,j-1,k) )/(2*dy);

                        // if(i==0){

                        //     phix[comp] = ( phi(comp,i+1,j,k) - phi(comp,nx-1,j,k) )/(2*dx);

                        // } else if(i==nx-1){

                        //     phix[comp] = ( phi(comp,0,j,k) - phi(comp,i-1,j,k) )/(2*dx);

                        // } else{

                        //     phix[comp] = ( phi(comp,i+1,j,k) - phi(comp,i-1,j,k) )/(2*dx);

                        // }

                        // if(j==0){

                        //     phiy[comp] = ( phi(comp,i,j+1,k) - phi(comp,i,ny-1,k) )/(2*dy);

                        // } else if(j==ny-1){

                        //     phiy[comp] = ( phi(comp,i,0,k) - phi(comp,i,j-1,k) )/(2*dy);

                        // } else{

                        //     phiy[comp] = ( phi(comp,i,j+1,k) - phi(comp,i,j-1,k) )/(2*dy);

                        // }

                        if(k==0){

                            phiz[comp] = ( phi(comp,i,j,k+1) - phi(comp,i,j,nz-1) )/(2*dz);

                        } else if(k==nz-1){

                            phiz[comp] = ( phi(comp,i,j,0) - phi(comp,i,j,k-1) )/(2*dz);

                        } else{

                            phiz[comp] = ( phi(comp,i,j,k+1) - phi(comp,i,j,k-1) )/(2*dz);

                        }

                    }

                    // Rolling back evolution of phit so timesteps align for energy calculation.

                    energydensity(i,j,k) = pow(phit(0,i,j,k) - 0.5*phitt(0,i,j,k)*dt,2) + pow(phit(1,i,j,k) - 0.5*phitt(1,i,j,k)*dt,2) 
                                         + pow(phix[0],2) + pow(phix[1],2) + pow(phiy[0],2) + pow(phiy[1],2) + pow(phiz[0],2) + pow(phiz[1],2) + 0.5*pow(phiMagSqr-1,2);

                    energy += dx*dy*dz*energydensity(i,j,k);


                    // Update field array

                    for(comp=0;comp<2;comp++){

                        phi(comp,i,j,k) += phit(comp,i,j,k)*dt;

                    }

                }
            }
        }


        valsPerLoop << energy << endl;


        if(makeGif && TimeStep%saveFreq == 0){

            ss.str(string());
            ss << gifFrame;
            string gifDataPath = dir_path + "/GifData/gifData_" + ss.str() + ".txt";
            ofstream gifData (gifDataPath.c_str());
            gifFrame+=1;

            for(i=0;i<nx;i++){
                for(j=0;j<ny;j++){
                    for(k=0;k<nz;k++){

                        gifData << pow(phi(0,i,j,k),2) + pow(phi(1,i,j,k),2) << endl;

                    }
                }
            }

        }

    }

    cout << endl;

    for(i=0;i<nx;i++){
        for(j=0;j<ny;j++){
            for(k=0;k<nz;k++){

                finalField << phi(0,i,j,k) << " " << phi(1,i,j,k) << endl;

            }
        }
    }


    gettimeofday(&end,NULL);

    cout << "Time taken: " << end.tv_sec - start.tv_sec << "s" << endl;


	return 0;

}


