#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <sys/time.h>
#include <omp.h>
#include <complex>

#include "header.hpp"

using namespace std;

// Perturb a straight, global string solution, evolve the system and investigate how it radiates.

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                 		  Parameters & Declarations
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

const int nx = 201;
const int ny = 201;
const int nz = 51;
const int nt = 1000;
const double dx = 0.7;
const double dy = 0.7;
const double dz = 0.7;
const double dt = 0.3;

const bool makeGif = false;
const int saveFreq = 10;


int main(){

	Array phi(2,2,nx,ny,nz,0.0), phitt(2,nx,ny,nz,0.0), bUpdate(2,nx,ny,nz,0.0), energydensity(nx,ny,nz,0.0);
	int comp, i, j, k, TimeStep, gifFrame, tNow, tPast, s;
    double phixx, phiyy, phizz, phiMagSqr, phix[2], phiy[2], phiz[2], phit[2], energy;

	struct timeval start, end;
    gettimeofday(&start, NULL);

    string file_path = __FILE__;
    string dir_path = file_path.substr(0,file_path.find_last_of('/'));
    stringstream ss;

    string icPath = dir_path + "/ic.txt";
    string finalFieldPath = dir_path + "/finalField.txt";
    string valsPerLoopPath = dir_path + "/valsPerLoop.txt";
    string testPath = dir_path + "/test.txt";

    ifstream ic (icPath.c_str());
    ofstream finalField (finalFieldPath.c_str());
    ofstream valsPerLoop (valsPerLoopPath.c_str());
    ofstream test (testPath.c_str());

    for(i=0;i<nx;i++){
        for(j=0;j<ny;j++){
            for(k=0;k<nz;k++){

                ic >> phi(0,0,i,j,k) >> phi(1,0,i,j,k);

                // Second time step is equal to first.

                phi(0,1,i,j,k) = phi(0,0,i,j,k);
                phi(1,1,i,j,k) = phi(1,0,i,j,k);

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

        tNow = (TimeStep+1)%2;
        tPast = TimeStep%2;

        // Set boundary conditions. Different loop sizes so that corner allocation is done correctly.

        // Sets Neumann boundary conditions on x and y as simple first attempt. Will need to use absorbing boundary conditions in future. z direction has periodic boundary conditions instead.

        for(i=1;i<nx-1;i++){
            for(k=1;k<nz-1;k++){

                phi(0,tNow,i,0,k) = phi(0,tNow,i,1,k);
                phi(1,tNow,i,0,k) = phi(1,tNow,i,1,k);

                phi(0,tNow,i,ny-1,k) = phi(0,tNow,i,ny-2,k);
                phi(1,tNow,i,ny-1,k) = phi(1,tNow,i,ny-2,k);

            }
        }

        for(j=0;j<ny;j++){
            for(k=1;k<nz-1;k++){

                phi(0,tNow,0,j,k) = phi(0,tNow,1,j,k);
                phi(1,tNow,0,j,k) = phi(1,tNow,1,j,k);

                phi(0,tNow,nx-1,j,k) = phi(0,tNow,nx-2,j,k);
                phi(1,tNow,nx-1,j,k) = phi(1,tNow,nx-2,j,k);

            }
        }


        // Calculate time derivatives using EoMs

        energy=0;

        #pragma omp parallel for reduction(+:energy) default(none) shared(phitt,phi,energydensity,tNow,tPast) private(phiMagSqr,phixx,phiyy,phizz,j,k,comp,phix,phiy,phiz,phit)
        for(i=1;i<nx-1;i++){
            for(j=1;j<ny-1;j++){
                for(k=1;k<nz-1;k++){

                    phiMagSqr = pow(phi(0,tNow,i,j,k),2) + pow(phi(1,tNow,i,j,k),2);

                    for(comp=0;comp<2;comp++){

                         // 2nd order spatial derivatives calculated with 2nd order finite difference

                        phixx = ( phi(comp,tNow,i+1,j,k) - 2*phi(comp,tNow,i,j,k) + phi(comp,tNow,i-1,j,k) )/(dx*dx);
                        phiyy = ( phi(comp,tNow,i,j+1,k) - 2*phi(comp,tNow,i,j,k) + phi(comp,tNow,i,j-1,k) )/(dy*dy);

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

                            phizz = ( phi(comp,tNow,i,j,k+1) - 2*phi(comp,tNow,i,j,k) + phi(comp,tNow,i,j,nz-1) )/(dz*dz);

                        } else if(k==nz-1){

                            phizz = ( phi(comp,tNow,i,j,0) - 2*phi(comp,tNow,i,j,k) + phi(comp,tNow,i,j,k-1) )/(dz*dz);

                        } else{
                        
                            phizz = ( phi(comp,tNow,i,j,k+1) - 2*phi(comp,tNow,i,j,k) + phi(comp,tNow,i,j,k-1) )/(dz*dz);

                        }

                        // Calculate second order time derivatives

                        phitt(comp,i,j,k) = phixx + phiyy + phizz - (phiMagSqr - 1)*phi(comp,tNow,i,j,k);


                        // Calculate first order derivatives for energy


                        phix[comp] = ( phi(comp,tNow,i+1,j,k) - phi(comp,tNow,i-1,j,k) )/(2*dx);
                        phiy[comp] = ( phi(comp,tNow,i,j+1,k) - phi(comp,tNow,i,j-1,k) )/(2*dy);


                        if(k==0){

                            phiz[comp] = ( phi(comp,tNow,i,j,k+1) - phi(comp,tNow,i,j,nz-1) )/(2*dz);

                        } else if(k==nz-1){

                            phiz[comp] = ( phi(comp,tNow,i,j,0) - phi(comp,tNow,i,j,k-1) )/(2*dz);

                        } else{

                            phiz[comp] = ( phi(comp,tNow,i,j,k+1) - phi(comp,tNow,i,j,k-1) )/(2*dz);

                        }

                        // Updated field value - previous field value to get time derivative with central differencing

                        phit[comp] = ( 2*phi(comp,tNow,i,j,k) - 2*phi(comp,tPast,i,j,k) + dt*dt*phitt(comp,i,j,k) )/(2*dt);

                    }

                    energydensity(i,j,k) = pow(phit[0],2) + pow(phit[1],2) + pow(phix[0],2) + pow(phix[1],2) + pow(phiy[0],2) + pow(phiy[1],2) + pow(phiz[0],2) + pow(phiz[1],2)
                                         + 0.5*pow(phiMagSqr-1,2);

                    energy += dx*dy*dz*energydensity(i,j,k);

                }
            }
        }

        // Set absorbing boundary conditions (for next time step)

       // Set along x boundary (so that absorbing condition is true at grid point neighbouring the boundary)
       

        // for(j=1;j<ny-1;j++){
        //     for(k=0;k<nz;k++){
        //         for(comp=0;comp<2;comp++){

        //             // left x boundary

        //             phiyy = ( phi(comp,tNow,1,j+1,k) - 2*phi(comp,tNow,1,j,k) + phi(comp,tNow,1,j-1,k) )/(dy*dy);

        //             if(k==0){

        //                 phizz = ( phi(comp,tNow,1,j,k+1) - 2*phi(comp,tNow,1,j,k) + phi(comp,tNow,1,j,nz-1) )/(dz*dz);

        //             } else if(k==nz-1){

        //                 phizz = ( phi(comp,tNow,1,j,0) - 2*phi(comp,tNow,1,j,k) + phi(comp,tNow,1,j,k-1) )/(dz*dz);

        //             } else{
                        
        //                 phizz = ( phi(comp,tNow,1,j,k+1) - 2*phi(comp,tNow,1,j,k) + phi(comp,tNow,1,j,k-1) )/(dz*dz);

        //             }

        //             bUpdate(comp,0,j,k) = 2*phi(comp,tNow,2,j,k) - 2*phi(comp,tPast,2,j,k) + dt*dt*phitt(comp,2,j,k) + phi(comp,tPast,0,j,k)
        //                                 - 4*dt*dx*( phitt(comp,1,j,k) - 0.5*(phiyy + phizz) );

        //             // right x boundary

        //             phiyy = ( phi(comp,tNow,nx-2,j+1,k) - 2*phi(comp,tNow,nx-2,j,k) + phi(comp,tNow,nx-2,j-1,k) )/(dy*dy);

        //             if(k==0){

        //                 phizz = ( phi(comp,tNow,nx-2,j,k+1) - 2*phi(comp,tNow,nx-2,j,k) + phi(comp,tNow,nx-2,j,nz-1) )/(dz*dz);

        //             } else if(k==nz-1){

        //                 phizz = ( phi(comp,tNow,nx-2,j,0) - 2*phi(comp,tNow,nx-2,j,k) + phi(comp,tNow,nx-2,j,k-1) )/(dz*dz);

        //             } else{
                        
        //                 phizz = ( phi(comp,tNow,nx-2,j,k+1) - 2*phi(comp,tNow,nx-2,j,k) + phi(comp,tNow,nx-2,j,k-1) )/(dz*dz);

        //             }

        //             bUpdate(comp,nx-1,j,k) = 2*phi(comp,tNow,nx-3,j,k) - 2*phi(comp,tPast,nx-3,j,k) + dt*dt*phitt(comp,nx-3,j,k) + phi(comp,tPast,nx-1,j,k)
        //                                    + 4*dt*dx*( phitt(comp,nx-2,j,k) - 0.5*(phiyy + phizz) );


            
        //         }
        //     }
        // } 

        // // Set along y boundary


        // for(i=1;i<nx-1;i++){
        //     for(k=0;k<nz;k++){
        //         for(comp=0;comp<2;comp++){

        //             // lower y boundary

        //             phixx = ( phi(comp,tNow,i+1,1,k) - 2*phi(comp,tNow,i,1,k) + phi(comp,tNow,i-1,1,k) )/(dx*dx);

        //             if(k==0){

        //                 phizz = ( phi(comp,tNow,i,1,k+1) - 2*phi(comp,tNow,i,1,k) + phi(comp,tNow,i,1,nz-1) )/(dz*dz);

        //             } else if(k==nz-1){

        //                 phizz = ( phi(comp,tNow,i,1,0) - 2*phi(comp,tNow,i,1,k) + phi(comp,tNow,i,1,k-1) )/(dz*dz);

        //             } else{
                        
        //                 phizz = ( phi(comp,tNow,i,1,k+1) - 2*phi(comp,tNow,i,1,k) + phi(comp,tNow,i,1,k-1) )/(dz*dz);

        //             }

        //             bUpdate(comp,i,0,k) = 2*phi(comp,tNow,i,2,k) - 2*phi(comp,tPast,i,2,k) + dt*dt*phitt(comp,i,2,k) + phi(comp,tPast,i,0,k)
        //                                 - 4*dt*dy*( phitt(comp,i,1,k) - 0.5*(phixx + phizz) );

        //             // upper y boundary

        //             phixx = ( phi(comp,tNow,i+1,ny-2,k) - 2*phi(comp,tNow,i,ny-2,k) + phi(comp,tNow,i-1,ny-2,k) )/(dx*dx);

        //             if(k==0){

        //                 phizz = ( phi(comp,tNow,i,ny-2,k+1) - 2*phi(comp,tNow,i,ny-2,k) + phi(comp,tNow,i,ny-2,nz-1) )/(dz*dz);

        //             } else if(k==nz-1){

        //                 phizz = ( phi(comp,tNow,i,ny-2,0) - 2*phi(comp,tNow,i,ny-2,k) + phi(comp,tNow,i,ny-2,k-1) )/(dz*dz);

        //             } else{
                        
        //                 phizz = ( phi(comp,tNow,i,ny-2,k+1) - 2*phi(comp,tNow,i,ny-2,k) + phi(comp,tNow,i,ny-2,k-1) )/(dz*dz);

        //             }

        //             bUpdate(comp,i,ny-1,k) = 2*phi(comp,tNow,i,ny-3,k) - 2*phi(comp,tPast,i,ny-3,k) + dt*dt*phitt(comp,i,ny-3,k) + phi(comp,tPast,i,ny-1,k)
        //                                    + 4*dt*dy*( phitt(comp,i,ny-2,k) - 0.5*(phixx + phizz) );

                
        //         }
        //     }

        // }


        // Absorbing boundary conditions so that equation is satisfied at the boundary.

        // for(j=1;j<ny-1;j++){
        //     for(k=0;k<nz;k++){
        //         for(comp=0;comp<2;comp++){

        //             // left x boundary

        //             phiyy = ( phi(comp,tNow,0,j+1,k) - 2*phi(comp,tNow,0,j,k) + phi(comp,tNow,0,j-1,k) )/(dy*dy);

        //             if(k==0){

        //                 phizz = ( phi(comp,tNow,0,j,k+1) - 2*phi(comp,tNow,0,j,k) + phi(comp,tNow,0,j,nz-1) )/(dz*dz);

        //             } else if(k==nz-1){

        //                 phizz = ( phi(comp,tNow,0,j,0) - 2*phi(comp,tNow,0,j,k) + phi(comp,tNow,0,j,k-1) )/(dz*dz);

        //             } else{
                        
        //                 phizz = ( phi(comp,tNow,0,j,k+1) - 2*phi(comp,tNow,0,j,k) + phi(comp,tNow,0,j,k-1) )/(dz*dz);

        //             }

        //             bUpdate(comp,0,j,k) = 2*phi(comp,tNow,2,j,k) - 2*phi(comp,tPast,2,j,k) + dt*dt*phitt(comp,2,j,k) + phi(comp,tPast,0,j,k)
        //                                 - 4*dt*dx*( phitt(comp,1,j,k) - 0.5*(phiyy + phizz) );

        //             // right x boundary

        //             phiyy = ( phi(comp,tNow,nx-2,j+1,k) - 2*phi(comp,tNow,nx-2,j,k) + phi(comp,tNow,nx-2,j-1,k) )/(dy*dy);

        //             if(k==0){

        //                 phizz = ( phi(comp,tNow,nx-2,j,k+1) - 2*phi(comp,tNow,nx-2,j,k) + phi(comp,tNow,nx-2,j,nz-1) )/(dz*dz);

        //             } else if(k==nz-1){

        //                 phizz = ( phi(comp,tNow,nx-2,j,0) - 2*phi(comp,tNow,nx-2,j,k) + phi(comp,tNow,nx-2,j,k-1) )/(dz*dz);

        //             } else{
                        
        //                 phizz = ( phi(comp,tNow,nx-2,j,k+1) - 2*phi(comp,tNow,nx-2,j,k) + phi(comp,tNow,nx-2,j,k-1) )/(dz*dz);

        //             }

        //             bUpdate(comp,nx-1,j,k) = 2*phi(comp,tNow,nx-3,j,k) - 2*phi(comp,tPast,nx-3,j,k) + dt*dt*phitt(comp,nx-3,j,k) + phi(comp,tPast,nx-1,j,k)
        //                                    + 4*dt*dx*( phitt(comp,nx-2,j,k) - 0.5*(phiyy + phizz) );


            
        //         }
        //     }
        // } 

        // // Set along y boundary


        // for(i=1;i<nx-1;i++){
        //     for(k=0;k<nz;k++){
        //         for(comp=0;comp<2;comp++){

        //             // lower y boundary

        //             phixx = ( phi(comp,tNow,i+1,1,k) - 2*phi(comp,tNow,i,1,k) + phi(comp,tNow,i-1,1,k) )/(dx*dx);

        //             if(k==0){

        //                 phizz = ( phi(comp,tNow,i,1,k+1) - 2*phi(comp,tNow,i,1,k) + phi(comp,tNow,i,1,nz-1) )/(dz*dz);

        //             } else if(k==nz-1){

        //                 phizz = ( phi(comp,tNow,i,1,0) - 2*phi(comp,tNow,i,1,k) + phi(comp,tNow,i,1,k-1) )/(dz*dz);

        //             } else{
                        
        //                 phizz = ( phi(comp,tNow,i,1,k+1) - 2*phi(comp,tNow,i,1,k) + phi(comp,tNow,i,1,k-1) )/(dz*dz);

        //             }

        //             bUpdate(comp,i,0,k) = 2*phi(comp,tNow,i,2,k) - 2*phi(comp,tPast,i,2,k) + dt*dt*phitt(comp,i,2,k) + phi(comp,tPast,i,0,k)
        //                                 - 4*dt*dy*( phitt(comp,i,1,k) - 0.5*(phixx + phizz) );

        //             // upper y boundary

        //             phixx = ( phi(comp,tNow,i+1,ny-2,k) - 2*phi(comp,tNow,i,ny-2,k) + phi(comp,tNow,i-1,ny-2,k) )/(dx*dx);

        //             if(k==0){

        //                 phizz = ( phi(comp,tNow,i,ny-2,k+1) - 2*phi(comp,tNow,i,ny-2,k) + phi(comp,tNow,i,ny-2,nz-1) )/(dz*dz);

        //             } else if(k==nz-1){

        //                 phizz = ( phi(comp,tNow,i,ny-2,0) - 2*phi(comp,tNow,i,ny-2,k) + phi(comp,tNow,i,ny-2,k-1) )/(dz*dz);

        //             } else{
                        
        //                 phizz = ( phi(comp,tNow,i,ny-2,k+1) - 2*phi(comp,tNow,i,ny-2,k) + phi(comp,tNow,i,ny-2,k-1) )/(dz*dz);

        //             }

        //             bUpdate(comp,i,ny-1,k) = 2*phi(comp,tNow,i,ny-3,k) - 2*phi(comp,tPast,i,ny-3,k) + dt*dt*phitt(comp,i,ny-3,k) + phi(comp,tPast,i,ny-1,k)
        //                                    + 4*dt*dy*( phitt(comp,i,ny-2,k) - 0.5*(phixx + phizz) );

                
        //         }
        //     }

        // }

        // // Now update the boundary values

        // for(i=0;i<nx;i=i+nx-1){
        //     for(j=1;j<ny-1;j++){
        //         for(k=0;k<nz;k++){
        //             for(comp=0;comp<2;comp++){

        //                 phi(comp,tPast,i,j,k) = bUpdate(comp,i,j,k);

        //             }
        //         }
        //     }
        // }

        // for(i=1;i<nx-1;i++){
        //     for(j=0;j<ny;j=j+ny-1){
        //         for(k=0;k<nz;k++){
        //             for(comp=0;comp<2;comp++){

        //                 phi(comp,tPast,i,j,k) = bUpdate(comp,i,j,k);

        //             }
        //         }
        //     }
        // }

        // May need to do this for the corner points too but I don't think they are actually ever used.


        // Update the array

        #pragma omp parallel for default(none) shared(phi,phitt,tNow,tPast) private(j,k,comp)
        for(i=1;i<nx-1;i++){
            for(j=1;j<ny-1;j++){
                for(k=0;k<nz;k++){
                    for(comp=0;comp<2;comp++){

                        phi(comp,tPast,i,j,k) = 2*phi(comp,tNow,i,j,k) - phi(comp,tPast,i,j,k) + dt*dt*phitt(comp,i,j,k);

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

                        gifData << pow(phi(0,tPast,i,j,k),2) + pow(phi(1,tPast,i,j,k),2) << endl;

                    }
                }
            }

        }

    }

    cout << "\rTimestep " << nt << " completed." << endl;

    for(i=0;i<nx;i++){
        for(j=0;j<ny;j++){
            for(k=0;k<nz;k++){

                finalField << phi(0,tPast,i,j,k) << " " << phi(1,tPast,i,j,k) << endl;

            }
        }
    }

    for(j=0;j<ny;j++){
        for(k=0;k<nz;k++){

            test << phi(0,tPast,0,j,k) << " " << phi(1,tPast,0,j,k) << endl;

        }
    }


    gettimeofday(&end,NULL);

    cout << "Time taken: " << end.tv_sec - start.tv_sec << "s" << endl;


	return 0;

}


