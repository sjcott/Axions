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
const int nt = 4000;
const double dx = 0.7;
const double dy = 0.7;
const double dz = 0.7;
const double dt = 0.3;

const bool makeGif = false;
const int saveFreq = 10;
const int countRate = 10;

const string xyBC = "absorbing";


int main(){

	Array phi(2,2,nx,ny,nz,0.0), phitt(2,nx,ny,nz,0.0), bUpdate(2,nx,ny,nz,0.0), energydensity(nx,ny,nz,0.0);
	int comp, i, j, k, TimeStep, gifFrame, tNow, tPast, s, counter;
    double phixx, phiyy, phizz, phiMagSqr, phix[2], phiy[2], phiz[2], phit[2], energy, phitx, phity;

	struct timeval start, end;
    gettimeofday(&start, NULL);

    string file_path = __FILE__;
    string dir_path = file_path.substr(0,file_path.find_last_of('/'));
    stringstream ss;

    string icPath = dir_path + "/ic.txt";
    string finalFieldPath = dir_path + "/finalField.txt";
    string valsPerLoopPath = dir_path + "/valsPerLoop.txt";
    string test1Path = dir_path + "/test1.txt";
    string test2Path = dir_path + "/test2.txt";

    ifstream ic (icPath.c_str());
    ofstream finalField (finalFieldPath.c_str());
    ofstream valsPerLoop (valsPerLoopPath.c_str());
    ofstream test1 (test1Path.c_str());
    ofstream test2 (test2Path.c_str());

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
    counter = 0;

    for(TimeStep=0;TimeStep<nt;TimeStep++){

        //if((1000*TimeStep)%(nt-1)==0){
        if(TimeStep>counter){

            cout << "\rTimestep " << TimeStep-1 << " completed." << flush;

            counter += countRate;

        }

        tNow = (TimeStep+1)%2;
        tPast = TimeStep%2;

        // Set boundary conditions. Different loop sizes so that corner allocation is done correctly.

        // Sets Neumann boundary conditions on x and y as simple first attempt. Will need to use absorbing boundary conditions in future. z direction has periodic boundary conditions instead.

        if(xyBC=="neumann"){

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

        //////////////////////////////////////////////////////////////////////////////////////////
        //                            Absorbing boundary conditions
        //////////////////////////////////////////////////////////////////////////////////////////

        if(xyBC=="absorbing"){

            // Set along x boundary

            // for(j=1;j<ny-1;j++){
            //     for(k=1;k<nz-1;k++){
            //         for(comp=0;comp<2;comp++){

            //             // Using points next to boundary to satisfy condition

            //             // x<0 boundary

            //             phiyy = ( phi(comp,tNow,1,j+1,k) - 2*phi(comp,tNow,1,j,k) + phi(comp,tNow,1,j-1,k) )/(dy*dy);
            //             phizz = ( phi(comp,tNow,1,j,k+1) - 2*phi(comp,tNow,1,j,k) + phi(comp,tNow,1,j,k-1) )/(dz*dz);

            //             bUpdate(comp,0,j,k) = phi(comp,tPast,0,j,k) + 2*phi(comp,tNow,2,j,k) - 2*phi(comp,tPast,2,j,k) + dt*dt*phitt(comp,2,j,k)
            //                                 - 4*dt*dx*( -phitt(comp,1,j,k) + 0.5*(phiyy + phizz) );

            //             // x>0 boundary

            //             phiyy = ( phi(comp,tNow,nx-2,j+1,k) - 2*phi(comp,tNow,nx-2,j,k) + phi(comp,tNow,nx-2,j-1,k) )/(dy*dy);
            //             phizz = ( phi(comp,tNow,nx-2,j,k+1) - 2*phi(comp,tNow,nx-2,j,k) + phi(comp,tNow,nx-2,j,k-1) )/(dz*dz);

            //             bUpdate(comp,nx-1,j,k) = phi(comp,tPast,nx-1,j,k) + 2*phi(comp,tNow,nx-3,j,k) - 2*phi(comp,tPast,nx-3,j,k) + dt*dt*phitt(comp,nx-3,j,k)
            //                                    - 4*dt*dx*( -phitt(comp,nx-2,j,k) + 0.5*(phiyy + phizz) );

            //         }
            //     }
            // }

            // // Set along y boundary

            // for(i=1;i<nx-1;i++){
            //     for(k=1;k<nz-1;k++){
            //         for(comp=0;comp<2;comp++){

            //             // Using points next to boundary again

            //             // y<0 boundary

            //             phixx = ( phi(comp,tNow,i+1,1,k) - 2*phi(comp,tNow,i,1,k) + phi(comp,tNow,i-1,1,k) )/(dx*dx);
            //             phizz = ( phi(comp,tNow,i,1,k+1) - 2*phi(comp,tNow,i,1,k) + phi(comp,tNow,i,1,k-1) )/(dz*dz);

            //             bUpdate(comp,i,0,k) = phi(comp,tPast,i,0,k) + 2*phi(comp,tNow,i,2,k) - 2*phi(comp,tPast,i,2,k) + dt*dt*phitt(comp,i,2,k)
            //                                 - 4*dt*dy*( -phitt(comp,i,1,k) + 0.5*(phixx + phizz) );

            //             // y>0 boundary

            //             phixx = ( phi(comp,tNow,i+1,ny-2,k) - 2*phi(comp,tNow,i,ny-2,k) + phi(comp,tNow,i-1,ny-2,k) )/(dx*dx);
            //             phizz = ( phi(comp,tNow,i,ny-2,k+1) - 2*phi(comp,tNow,i,ny-2,k) + phi(comp,tNow,i,ny-2,k-1) )/(dz*dz);

            //             bUpdate(comp,i,ny-1,k) = phi(comp,tPast,i,ny-1,k) + 2*phi(comp,tNow,i,ny-3,k) - 2*phi(comp,tPast,i,ny-3,k) + dt*dt*phitt(comp,i,ny-3,k)
            //                                    - 4*dt*dy*( -phitt(comp,i,ny-2,k) + 0.5*(phixx + phizz) );

            //         }
            //     }
            // }


            // // Now need to update the boundary values. Corners have been ignored because they don't have any impact on dynamics with this method

            // for(j=1;j<ny-1;j++){
            //     for(k=1;k<nz-1;k++){
            //         for(comp=0;comp<2;comp++){

            //             phi(comp,tPast,0,j,k) = bUpdate(comp,0,j,k);
            //             phi(comp,tPast,nx-1,j,k) = bUpdate(comp,nx-1,j,k);

            //         }
            //     }
            // }

            // for(i=1;i<nx-1;i++){
            //     for(k=1;k<nz-1;k++){
            //         for(comp=0;comp<2;comp++){

            //             phi(comp,tPast,i,0,k) = bUpdate(comp,i,0,k);
            //             phi(comp,tPast,i,ny-1,k) = bUpdate(comp,i,ny-1,k);

            //         }
            //     }
            // }

            // x boundary

            for(j=1;j<ny-1;j++){
                for(k=1;k<nz-1;k++){
                    for(comp=0;comp<2;comp++){

                        // x<0 boundary

                        phiyy = ( phi(comp,tNow,0,j+1,k) - 2*phi(comp,tNow,0,j,k) + phi(comp,tNow,0,j-1,k) )/(dy*dy);
                        phizz = ( phi(comp,tNow,0,j,k+1) - 2*phi(comp,tNow,0,j,k) + phi(comp,tNow,0,j,k-1) )/(dz*dz);

                        phitx = ( phi(comp,tNow,1,j,k) - phi(comp,tNow,0,j,k) - phi(comp,tPast,1,j,k) + phi(comp,tPast,0,j,k) )/(dt*dx);

                        phitt(comp,0,j,k) = phitx + 0.5*(phiyy + phizz);

                        // x>0 boundary

                        phiyy = ( phi(comp,tNow,nx-1,j+1,k) - 2*phi(comp,tNow,nx-1,j,k) + phi(comp,tNow,nx-1,j-1,k) )/(dy*dy);
                        phizz = ( phi(comp,tNow,nx-1,j,k+1) - 2*phi(comp,tNow,nx-1,j,k) + phi(comp,tNow,nx-1,j,k-1) )/(dz*dz);

                        phitx = ( phi(comp,tNow,nx-1,j,k) - phi(comp,tNow,nx-2,j,k) - phi(comp,tPast,nx-1,j,k) + phi(comp,tPast,nx-2,j,k) )/(dt*dx);

                        phitt(comp,nx-1,j,k) = -phitx + 0.5*(phiyy + phizz);

                    }

                    if(TimeStep==nt-1){test1 << phitt(0,0,j,k) << " " << phitt(1,0,j,k) << " " << phitt(0,nx-1,j,k) << " " << phitt(1,nx-1,j,k) << endl;}

                }
            }

            // Set along y boundary

            for(i=1;i<nx-1;i++){
                for(k=1;k<nz-1;k++){
                    for(comp=0;comp<2;comp++){

                        // y<0 boundary

                        phixx = ( phi(comp,tNow,i+1,0,k) - 2*phi(comp,tNow,i,0,k) + phi(comp,tNow,i-1,0,k) )/(dx*dx);
                        phizz = ( phi(comp,tNow,i,0,k+1) - 2*phi(comp,tNow,i,0,k) + phi(comp,tNow,i,0,k-1) )/(dz*dz);

                        phity = ( phi(comp,tNow,i,1,k) - phi(comp,tNow,i,0,k) - phi(comp,tPast,i,1,k) + phi(comp,tPast,i,0,k) )/(dt*dy);

                        phitt(comp,i,0,k) = phity + 0.5*(phixx + phizz);


                        // y>0 boundary

                        phixx = ( phi(comp,tNow,i+1,ny-1,k) - 2*phi(comp,tNow,i,ny-1,k) + phi(comp,tNow,i-1,ny-1,k) )/(dx*dx);
                        phizz = ( phi(comp,tNow,i,ny-1,k+1) - 2*phi(comp,tNow,i,ny-1,k) + phi(comp,tNow,i,ny-1,k-1) )/(dz*dz);

                        phity = ( phi(comp,tNow,i,ny-1,k) - phi(comp,tNow,i,ny-2,k) - phi(comp,tPast,i,ny-1,k) + phi(comp,tPast,i,ny-2,k) )/(dt*dy);

                        phitt(comp,i,ny-1,k) = -phity + 0.5*(phixx + phizz);

                    }

                    if(TimeStep==nt-1){test2 << phitt(0,i,0,k) << " " << phitt(1,i,0,k) << " " << phitt(0,i,ny-1,k) << " " << phitt(1,i,ny-1,k) << endl;}

                }
            }

            // Assign corners. Calculated by adding both absorbing both condition equations and subtracting 1/2 times the wave equation

            for(k=1;k<nz-1;k++){
                for(comp=0;comp<2;comp++){

                    // x,y<0 corner

                    phizz = ( phi(comp,tNow,0,0,k+1) - 2*phi(comp,tNow,0,0,k) + phi(comp,tNow,0,0,k-1) )/(dz*dz);

                    phitx = ( phi(comp,tNow,1,0,k) - phi(comp,tNow,0,0,k) - phi(comp,tPast,1,0,k) + phi(comp,tPast,0,0,k) )/(dt*dx);
                    phity = ( phi(comp,tNow,0,1,k) - phi(comp,tNow,0,0,k) - phi(comp,tPast,0,1,k) + phi(comp,tPast,0,0,k) )/(dt*dy);

                    phitt(comp,0,0,k) = ( phizz + 2*(phitx + phity) )/3;

                    // x<0,y>0 corner

                    phizz = ( phi(comp,tNow,0,ny-1,k+1) - 2*phi(comp,tNow,0,ny-1,k) + phi(comp,tNow,0,ny-1,k-1) )/(dz*dz);

                    phitx = ( phi(comp,tNow,1,ny-1,k) - phi(comp,tNow,0,ny-1,k) - phi(comp,tPast,1,ny-1,k) + phi(comp,tPast,0,ny-1,k) )/(dt*dx);
                    phity = ( phi(comp,tNow,0,ny-1,k) - phi(comp,tNow,0,ny-2,k) - phi(comp,tPast,0,ny-1,k) + phi(comp,tPast,0,ny-2,k) )/(dt*dy);

                    phitt(comp,0,ny-1,k) = ( phizz + 2*(phitx - phity) )/3;

                    // x>0,y<0 corner

                    phizz = ( phi(comp,tNow,nx-1,0,k+1) - 2*phi(comp,tNow,nx-1,0,k) + phi(comp,tNow,nx-1,0,k-1) )/(dz*dz);

                    phitx = ( phi(comp,tNow,nx-1,0,k) - phi(comp,tNow,nx-2,0,k) - phi(comp,tPast,nx-1,0,k) + phi(comp,tPast,nx-2,0,k) )/(dt*dx);
                    phity = ( phi(comp,tNow,nx-1,1,k) - phi(comp,tNow,nx-1,0,k) - phi(comp,tPast,nx-1,1,k) + phi(comp,tPast,nx-1,0,k) )/(dt*dy);

                    phitt(comp,nx-1,0,k) = ( phizz - 2*(phitx - phity) )/3;

                    // x,y>0 corner

                    phizz = ( phi(comp,tNow,nx-1,ny-1,k+1) - 2*phi(comp,tNow,nx-1,ny-1,k) + phi(comp,tNow,nx-1,ny-1,k-1) )/(dz*dz);

                    phitx = ( phi(comp,tNow,nx-1,ny-1,k) - phi(comp,tNow,nx-2,ny-1,k) - phi(comp,tPast,nx-1,ny-1,k) + phi(comp,tPast,nx-2,ny-1,k) )/(dt*dx);
                    phity = ( phi(comp,tNow,nx-1,ny-1,k) - phi(comp,tNow,nx-1,ny-2,k) - phi(comp,tPast,nx-1,ny-1,k) + phi(comp,tPast,nx-1,ny-2,k) )/(dt*dy);

                    phitt(comp,nx-1,ny-1,k) = ( phizz - 2*(phitx + phity) )/3;

                }
            }



            // Set along x boundary

            // for(j=1;j<ny-1;j++){
            //     for(k=1;k<nz-1;k++){
            //         for(comp=0;comp<2;comp++){

            //             // x<0 boundary

            //             phiyy = ( phi(comp,tNow,0,j+1,k) - 2*phi(comp,tNow,0,j,k) + phi(comp,tNow,0,j-1,k) )/(dy*dy);
            //             phizz = ( phi(comp,tNow,0,j,k+1) - 2*phi(comp,tNow,0,j,k) + phi(comp,tNow,0,j,k-1) )/(dz*dz);

            //             bUpdate(comp,0,j,k) = ( dt*( 2*phi(comp,tNow,1,j,k) - 2*phi(comp,tPast,1,j,k) + dt*dt*phitt(comp,1,j,k) + phi(comp,tPast,0,j,k) )
            //                                 + 2*dx*( 2*phi(comp,tNow,0,j,k) - phi(comp,tPast,0,j,k) ) + dt*dt*dx*( phiyy + phizz ) )/(dt + 2*dx);

            //             // x>0 boundary

            //             phiyy = ( phi(comp,tNow,nx-1,j+1,k) - 2*phi(comp,tNow,nx-1,j,k) + phi(comp,tNow,nx-1,j-1,k) )/(dy*dy);
            //             phizz = ( phi(comp,tNow,nx-1,j,k+1) - 2*phi(comp,tNow,nx-1,j,k) + phi(comp,tNow,nx-1,j,k-1) )/(dz*dz);

            //             bUpdate(comp,nx-1,j,k) = ( dt*( 2*phi(comp,tNow,nx-2,j,k) - 2*phi(comp,tPast,nx-2,j,k) + dt*dt*phitt(comp,nx-2,j,k) + phi(comp,tPast,nx-1,j,k) )
            //                                    - 2*dx*( 2*phi(comp,tNow,nx-1,j,k) - phi(comp,tPast,nx-1,j,k) ) - dt*dt*dx*( phiyy + phizz ) )/(dt - 2*dx);

            //         }
            //     }
            // }

            // // Set along y boundary with conditions for corners

            // for(i=0;i<nx;i++){
            //     for(k=1;k<nz-1;k++){
            //         for(comp=0;comp<2;comp++){

            //             // y < 0 boundary

            //             phizz = ( phi(comp,tNow,i,0,k+1) - 2*phi(comp,tNow,i,0,k) + phi(comp,tNow,i,0,k-1) )/(dz*dz);

            //             if(i==0){

            //                 phixx = ( phi(comp,tNow,i+2,0,k) - 2*phi(comp,tNow,i+1,0,k) + phi(comp,tNow,i,0,k) )/(dx*dx);

            //                 bUpdate(comp,i,0,k) = ( dt*( bUpdate(comp,i,1,k) - phi(comp,tPast,i,1,k) + phi(comp,tPast,i,0,k) ) + 2*dx*( 2*phi(comp,tNow,i,0,k) - phi(comp,tPast,i,0,k) )
            //                                     + dt*dt*dx*( phixx + phizz ) )/(dt + 2*dx);

            //             } else if(i==nx-1){

            //                 phixx = ( phi(comp,tNow,i,0,k) - 2*phi(comp,tNow,i-1,0,k) + phi(comp,tNow,i-2,0,k) )/(dx*dx);

            //                 // Same update as above but different for non-boundary terms
            //                 bUpdate(comp,i,0,k) = ( dt*( bUpdate(comp,i,1,k) - phi(comp,tPast,i,1,k) + phi(comp,tPast,i,0,k) ) + 2*dx*( 2*phi(comp,tNow,i,0,k) - phi(comp,tPast,i,0,k) )
            //                                     + dt*dt*dx*( phixx + phizz ) )/(dt + 2*dx);

            //             } else{

            //                 phixx = ( phi(comp,tNow,i+1,0,k) - 2*phi(comp,tNow,i,0,k) + phi(comp,tNow,i-1,0,k) )/(dx*dx);

            //                 bUpdate(comp,i,0,k) = ( dt*( 2*phi(comp,tNow,i,1,k) - 2*phi(comp,tPast,i,1,k) + dt*dt*phitt(comp,i,1,k) + phi(comp,tPast,i,0,k) )
            //                                     + 2*dx*( 2*phi(comp,tNow,i,0,k) - phi(comp,tPast,i,0,k) ) + dt*dt*dx*( phixx + phizz ) )/(dt + 2*dx);

            //             }

            //             // y > 0 boundary

            //             phizz = ( phi(comp,tNow,i,ny-1,k+1) - 2*phi(comp,tNow,i,ny-1,k) + phi(comp,tNow,i,ny-1,k-1) )/(dz*dz);

            //             if(i==0){

            //                 phixx = ( phi(comp,tNow,i+2,ny-1,k) - 2*phi(comp,tNow,i+1,ny-1,k) + phi(comp,tNow,i,ny-1,k) )/(dx*dx);

            //                 bUpdate(comp,i,ny-1,k) = ( dt*( bUpdate(comp,i,ny-2,k) - phi(comp,tPast,i,ny-2,k) + phi(comp,tPast,i,ny-1,k) ) - 2*dx*( 2*phi(comp,tNow,i,ny-1,k) - phi(comp,tPast,i,ny-1,k) )
            //                                        - dt*dt*dx*( phixx + phizz ) )/(dt - 2*dx);

            //             } else if(i==nx-1){

            //                 phixx = ( phi(comp,tNow,i,ny-1,k) - 2*phi(comp,tNow,i-1,ny-1,k) + phi(comp,tNow,i-2,ny-1,k) )/(dx*dx);

            //                 bUpdate(comp,i,ny-1,k) = ( dt*( bUpdate(comp,i,ny-2,k) - phi(comp,tPast,i,ny-2,k) + phi(comp,tPast,i,ny-1,k) ) - 2*dx*( 2*phi(comp,tNow,i,ny-1,k) - phi(comp,tPast,i,ny-1,k) )
            //                                        - dt*dt*dx*( phixx + phizz ) )/(dt - 2*dx);

            //             } else{

            //                 phixx = ( phi(comp,tNow,i+1,ny-1,k) - 2*phi(comp,tNow,i,ny-1,k) + phi(comp,tNow,i-1,ny-1,k) )/(dx*dx);

            //                 bUpdate(comp,i,ny-1,k) = ( dt*( 2*phi(comp,tNow,i,ny-2,k) - 2*phi(comp,tPast,i,ny-2,k) + dt*dt*phitt(comp,i,ny-2,k) + phi(comp,tPast,i,ny-1,k) )
            //                                        - 2*dx*( 2*phi(comp,tNow,i,ny-1,k) - phi(comp,tPast,i,ny-1,k) ) - dt*dt*dx*( phixx + phizz ) )/(dt - 2*dx);

            //             }

            //         }
            //     }
            // }



            // // Now need to update the boundary values

            // for(j=1;j<ny-1;j++){
            //     for(k=1;k<nz-1;k++){
            //         for(comp=0;comp<2;comp++){

            //             phi(comp,tPast,0,j,k) = bUpdate(comp,0,j,k);
            //             phi(comp,tPast,nx-1,j,k) = bUpdate(comp,nx-1,j,k);

            //         }
            //     }
            // }

            // for(i=0;i<nx;i++){
            //     for(k=1;k<nz-1;k++){
            //         for(comp=0;comp<2;comp++){

            //             phi(comp,tPast,i,0,k) = bUpdate(comp,i,0,k);
            //             phi(comp,tPast,i,ny-1,k) = bUpdate(comp,i,ny-1,k);

            //         }
            //     }
            // }

        }



        // Update the array

        // #pragma omp parallel for default(none) shared(phi,phitt,tNow,tPast) private(j,k,comp)
        // for(i=1;i<nx-1;i++){
        //     for(j=1;j<ny-1;j++){
        //         for(k=1;k<nz-1;k++){
        //             for(comp=0;comp<2;comp++){

        //                 phi(comp,tPast,i,j,k) = 2*phi(comp,tNow,i,j,k) - phi(comp,tPast,i,j,k) + dt*dt*phitt(comp,i,j,k);

        //             }
        //         }
        //     }
        // }

        if(xyBC == "absorbing"){

            #pragma omp parallel for default(none) shared(phi,phitt,tNow,tPast) private(j,k,comp)
            for(i=0;i<nx;i++){
                for(j=0;j<ny;j++){
                    for(k=1;k<nz-1;k++){
                        for(comp=0;comp<2;comp++){

                            phi(comp,tPast,i,j,k) = 2*phi(comp,tNow,i,j,k) - phi(comp,tPast,i,j,k) + dt*dt*phitt(comp,i,j,k);

                        }
                    }
                }
            }

        } else{

            #pragma omp parallel for default(none) shared(phi,phitt,tNow,tPast) private(j,k,comp)
            for(i=1;i<nx-1;i++){
                for(j=1;j<ny-1;j++){
                    for(k=1;k<nz-1;k++){
                        for(comp=0;comp<2;comp++){

                            phi(comp,tPast,i,j,k) = 2*phi(comp,tNow,i,j,k) - phi(comp,tPast,i,j,k) + dt*dt*phitt(comp,i,j,k);

                        }
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


    gettimeofday(&end,NULL);

    cout << "Time taken: " << end.tv_sec - start.tv_sec << "s" << endl;


	return 0;

}


