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
const int nt = 6001;
const double dx = 0.7;
const double dy = 0.7;
const double dz = 0.7;
const double dt = 0.3;

const double eta1 = 1;
const double eta2 = 1;
const double lambda1 = 2;
const double lambda2 = 2;
const double lambda12 = 0;

const int q1 = 2;
const int q2 = q1 - 1;
const double g = sqrt(0.5*lambda1*pow(eta1,2)/(q1*q1 + q2*q2)); // Specific choice that sets all masses equal (assuming m1=m2)

const bool makeGif = true;
const int saveFreq = 10;
const int countRate = 10;

const string xyBC = "fixed"; // Don't try to use absorbing boundary conditions as not rewritten code for this yet


int main(){

	Array phi1(2,2,nx,ny,nz,0.0), phi2(2,2,nx,ny,nz,0.0), theta(3,2,nx,ny,nz,0.0), phi1tt(2,nx,ny,nz,0.0), phi2tt(2,nx,ny,nz,0.0), thetatt(3,nx,ny,nz,0.0), energydensity(nx,ny,nz,0.0), gaussDeviation(nx,ny,nz,0.0);
	int comp, i, j, k, TimeStep, gifFrame, tNow, tPast, s, counter;
    double phi1xx, phi1yy, phi1zz, phi1MagSqr, phi1x[2], phi1y[2], phi1z[2], cur1x, cur1y, cur1z, Fxy_y, Fxz_z, Fyx_x, Fyz_z, Fzx_x, Fzy_y, phi1t[2], energy, phi1tx, phi1ty, thetat[3],
           divTheta[2], divThetat, phi2xx, phi2yy, phi2zz, phi2MagSqr, phi2x[2], phi2y[2], phi2z[2], cur2x, cur2y, cur2z, phi2t[2], phi2tx, phi2ty,
           thetaDotCont, Fxy, Fxz, Fyz, FCont, deviationParameter, thetaxx, thetayy, thetazz, thetatx, thetaty;
    int c[2] = {1,-1}; // Useful definition to allow covariant deviative to be calculated when looping over components.

	struct timeval start, end;
    gettimeofday(&start, NULL);

    string file_path = __FILE__;
    string dir_path = file_path.substr(0,file_path.find_last_of('/'));
    stringstream ss;

    string icPath = dir_path + "/Data/moore_ic.txt";
    string finalFieldPath = dir_path + "/Data/finalField.txt";
    string valsPerLoopPath = dir_path + "/Data/valsPerLoop.txt";
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

                ic >> phi1(0,0,i,j,k) >> phi1(1,0,i,j,k) >> phi2(0,0,i,j,k) >> phi2(1,0,i,j,k) >> theta(0,0,i,j,k) >> theta(1,0,i,j,k) >> theta(2,0,i,j,k);

                // Second time step is equal to first.

                phi1(0,1,i,j,k) = phi1(0,0,i,j,k);
                phi1(1,1,i,j,k) = phi1(1,0,i,j,k);
                phi2(0,1,i,j,k) = phi2(0,0,i,j,k);
                phi2(1,1,i,j,k) = phi2(1,0,i,j,k);
                theta(0,1,i,j,k) = theta(0,0,i,j,k);
                theta(1,1,i,j,k) = theta(1,0,i,j,k);
                theta(2,1,i,j,k) = theta(2,0,i,j,k);

            }
        }
    }

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

                    // Doing theta first because it is used in one of the covariant derivative of phi
                    for(comp=0;comp<3;comp++){

                        theta(comp,tNow,i,0,k) = theta(comp,tNow,i,1,k);
                        theta(comp,tNow,i,ny-1,k) = theta(comp,tNow,i,ny-2,k);

                    }

                    for(comp=0;comp<2;comp++){

                        phi1(comp,tNow,i,0,k) = cos(q1*theta(1,tNow,i,0,k))*phi1(comp,tNow,i,1,k) + c[comp]*sin(q1*theta(1,tNow,i,0,k))*phi1(comp+c[comp],tNow,i,1,k);
                        phi1(comp,tNow,i,ny-1,k) = cos(q1*theta(1,tNow,i,ny-2,k))*phi1(comp,tNow,i,ny-2,k) - c[comp]*sin(q1*theta(1,tNow,i,ny-2,k))*phi1(comp+c[comp],tNow,i,ny-2,k);

                        phi2(comp,tNow,i,0,k) = cos(q2*theta(1,tNow,i,0,k))*phi2(comp,tNow,i,1,k) + c[comp]*sin(q2*theta(1,tNow,i,0,k))*phi2(comp+c[comp],tNow,i,1,k);
                        phi2(comp,tNow,i,ny-1,k) = cos(q2*theta(1,tNow,i,ny-2,k))*phi2(comp,tNow,i,ny-2,k) - c[comp]*sin(q2*theta(1,tNow,i,ny-2,k))*phi2(comp+c[comp],tNow,i,ny-2,k);

                    }

                }
            }

            for(j=0;j<ny;j++){
                for(k=1;k<nz-1;k++){

                    // Doing theta first because it is used in one of the covariant derivative of phi
                    for(comp=0;comp<3;comp++){

                        theta(comp,tNow,0,j,k) = theta(comp,tNow,1,j,k);
                        theta(comp,tNow,nx-1,j,k) = theta(comp,tNow,nx-2,j,k);

                    }

                    for(comp=0;comp<2;comp++){

                        phi1(comp,tNow,0,j,k) = cos(q1*theta(0,tNow,0,j,k))*phi1(comp,tNow,1,j,k) + c[comp]*sin(q1*theta(0,tNow,0,j,k))*phi1(comp+c[comp],tNow,1,j,k);
                        phi1(comp,tNow,nx-1,j,k) = cos(q1*theta(0,tNow,nx-2,j,k))*phi1(comp,tNow,nx-2,j,k) - c[comp]*sin(q1*theta(0,tNow,nx-2,j,k))*phi1(comp+c[comp],tNow,nx-2,j,k);

                        phi2(comp,tNow,0,j,k) = cos(q2*theta(0,tNow,0,j,k))*phi2(comp,tNow,1,j,k) + c[comp]*sin(q2*theta(0,tNow,0,j,k))*phi2(comp+c[comp],tNow,1,j,k);
                        phi2(comp,tNow,nx-1,j,k) = cos(q2*theta(0,tNow,nx-2,j,k))*phi2(comp,tNow,nx-2,j,k) - c[comp]*sin(q2*theta(0,tNow,nx-2,j,k))*phi2(comp+c[comp],tNow,nx-2,j,k);

                    }

                    // phi(0,tNow,0,j,k) = phi(0,tNow,1,j,k);
                    // phi(1,tNow,0,j,k) = phi(1,tNow,1,j,k);
                    // theta(0,tNow,0,j,k) = theta(0,tNow,1,j,k);
                    // theta(1,tNow,0,j,k) = theta(1,tNow,1,j,k);
                    // theta(2,tNow,0,j,k) = theta(2,tNow,1,j,k);

                    // phi(0,tNow,nx-1,j,k) = phi(0,tNow,nx-2,j,k);
                    // phi(1,tNow,nx-1,j,k) = phi(1,tNow,nx-2,j,k);
                    // theta(0,tNow,nx-1,j,k) = theta(0,tNow,nx-2,j,k);
                    // theta(1,tNow,nx-1,j,k) = theta(1,tNow,nx-2,j,k);
                    // theta(2,tNow,nx-1,j,k) = theta(2,tNow,nx-2,j,k);

                }
            }

        }


        // Calculate time derivatives using EoMs

        energy=0;
        deviationParameter = 0;

        #pragma omp parallel for reduction(+:energy,deviationParameter) default(none) shared(phi1,phi2,theta,phi1tt,phi2tt,thetatt,energydensity,gaussDeviation,tNow,tPast,c,TimeStep,test2) \
                private(phi1MagSqr,phi2MagSqr,phi1xx,phi1yy,phi1zz,phi2xx,phi2yy,phi2zz,j,k,comp,phi1x,phi1y,phi1z,phi1t,cur1x,cur1y,cur1z,phi2x,phi2y,phi2z,phi2t,cur2x,cur2y,cur2z, \
                    Fxy_y,Fxz_z,Fyx_x,Fyz_z,Fzx_x,Fzy_y,divTheta,divThetat,thetat,thetaDotCont,Fxy,Fxz,Fyz,FCont)
        for(i=1;i<nx-1;i++){
            for(j=1;j<ny-1;j++){
                for(k=1;k<nz-1;k++){

                    phi1MagSqr = pow(phi1(0,tNow,i,j,k),2) + pow(phi1(1,tNow,i,j,k),2);
                    phi2MagSqr = pow(phi2(0,tNow,i,j,k),2) + pow(phi2(1,tNow,i,j,k),2);

                    // Loop over phi components

                    for(comp=0;comp<2;comp++){

                         // 2nd order spatial derivatives calculated with 2nd order finite difference

                        // c is 1 when comp = 0 and -1 when comp = 1

                        phi1xx = ( cos(q1*theta(0,tNow,i,j,k))*phi1(comp,tNow,i+1,j,k) + c[comp]*sin(q1*theta(0,tNow,i,j,k))*phi1(comp+c[comp],tNow,i+1,j,k) - 2*phi1(comp,tNow,i,j,k)
                                 + cos(q1*theta(0,tNow,i-1,j,k))*phi1(comp,tNow,i-1,j,k) - c[comp]*sin(q1*theta(0,tNow,i-1,j,k))*phi1(comp+c[comp],tNow,i-1,j,k) )/(dx*dx);

                        phi2xx = ( cos(q2*theta(0,tNow,i,j,k))*phi2(comp,tNow,i+1,j,k) + c[comp]*sin(q2*theta(0,tNow,i,j,k))*phi2(comp+c[comp],tNow,i+1,j,k) - 2*phi2(comp,tNow,i,j,k)
                                 + cos(q2*theta(0,tNow,i-1,j,k))*phi2(comp,tNow,i-1,j,k) - c[comp]*sin(q2*theta(0,tNow,i-1,j,k))*phi2(comp+c[comp],tNow,i-1,j,k) )/(dx*dx);


                        phi1yy = ( cos(q1*theta(1,tNow,i,j,k))*phi1(comp,tNow,i,j+1,k) + c[comp]*sin(q1*theta(1,tNow,i,j,k))*phi1(comp+c[comp],tNow,i,j+1,k) - 2*phi1(comp,tNow,i,j,k)
                                 + cos(q1*theta(1,tNow,i,j-1,k))*phi1(comp,tNow,i,j-1,k) - c[comp]*sin(q1*theta(1,tNow,i,j-1,k))*phi1(comp+c[comp],tNow,i,j-1,k) )/(dy*dy);

                        phi2yy = ( cos(q2*theta(1,tNow,i,j,k))*phi2(comp,tNow,i,j+1,k) + c[comp]*sin(q2*theta(1,tNow,i,j,k))*phi2(comp+c[comp],tNow,i,j+1,k) - 2*phi2(comp,tNow,i,j,k)
                                 + cos(q2*theta(1,tNow,i,j-1,k))*phi2(comp,tNow,i,j-1,k) - c[comp]*sin(q2*theta(1,tNow,i,j-1,k))*phi2(comp+c[comp],tNow,i,j-1,k) )/(dy*dy);


                        // enforce periodic boundary condition in z direction

                        if(k==0){

                            phi1zz = ( cos(q1*theta(2,tNow,i,j,k))*phi1(comp,tNow,i,j,k+1) + c[comp]*sin(q1*theta(2,tNow,i,j,k))*phi1(comp+c[comp],tNow,i,j,k+1) - 2*phi1(comp,tNow,i,j,k)
                                     + cos(q1*theta(2,tNow,i,j,nz-1))*phi1(comp,tNow,i,j,nz-1) - c[comp]*sin(q1*theta(2,tNow,i,j,nz-1))*phi1(comp+c[comp],tNow,i,j,nz-1) )/(dz*dz);

                            phi2zz = ( cos(q2*theta(2,tNow,i,j,k))*phi2(comp,tNow,i,j,k+1) + c[comp]*sin(q2*theta(2,tNow,i,j,k))*phi2(comp+c[comp],tNow,i,j,k+1) - 2*phi2(comp,tNow,i,j,k)
                                     + cos(q2*theta(2,tNow,i,j,nz-1))*phi2(comp,tNow,i,j,nz-1) - c[comp]*sin(q2*theta(2,tNow,i,j,nz-1))*phi2(comp+c[comp],tNow,i,j,nz-1) )/(dz*dz);

                        } else if(k==nz-1){

                            phi1zz = ( cos(q1*theta(2,tNow,i,j,k))*phi1(comp,tNow,i,j,0) + c[comp]*sin(q1*theta(2,tNow,i,j,k))*phi1(comp+c[comp],tNow,i,j,0) - 2*phi1(comp,tNow,i,j,k) 
                                     + cos(q1*theta(2,tNow,i,j,k-1))*phi1(comp,tNow,i,j,k-1) - c[comp]*sin(q1*theta(2,tNow,i,j,k-1))*phi1(comp+c[comp],tNow,i,j,k-1) )/(dz*dz);

                            phi2zz = ( cos(q2*theta(2,tNow,i,j,k))*phi2(comp,tNow,i,j,0) + c[comp]*sin(q2*theta(2,tNow,i,j,k))*phi2(comp+c[comp],tNow,i,j,0) - 2*phi2(comp,tNow,i,j,k) 
                                     + cos(q2*theta(2,tNow,i,j,k-1))*phi2(comp,tNow,i,j,k-1) - c[comp]*sin(q2*theta(2,tNow,i,j,k-1))*phi2(comp+c[comp],tNow,i,j,k-1) )/(dz*dz);

                        } else{
                        
                            phi1zz = ( cos(q1*theta(2,tNow,i,j,k))*phi1(comp,tNow,i,j,k+1) + c[comp]*sin(q1*theta(2,tNow,i,j,k))*phi1(comp+c[comp],tNow,i,j,k+1) - 2*phi1(comp,tNow,i,j,k) 
                                     + cos(q1*theta(2,tNow,i,j,k-1))*phi1(comp,tNow,i,j,k-1) - c[comp]*sin(q1*theta(2,tNow,i,j,k-1))*phi1(comp+c[comp],tNow,i,j,k-1) )/(dz*dz);

                            phi2zz = ( cos(q2*theta(2,tNow,i,j,k))*phi2(comp,tNow,i,j,k+1) + c[comp]*sin(q2*theta(2,tNow,i,j,k))*phi2(comp+c[comp],tNow,i,j,k+1) - 2*phi2(comp,tNow,i,j,k) 
                                     + cos(q2*theta(2,tNow,i,j,k-1))*phi2(comp,tNow,i,j,k-1) - c[comp]*sin(q2*theta(2,tNow,i,j,k-1))*phi2(comp+c[comp],tNow,i,j,k-1) )/(dz*dz);

                        }

                        // Calculate second order time derivatives

                        phi1tt(comp,i,j,k) = phi1xx + phi1yy + phi1zz - ( 0.5*lambda1*(phi1MagSqr - pow(eta1,2)) + 0.5*lambda12*(phi2MagSqr - pow(eta2,2)) )*phi1(comp,tNow,i,j,k);
                        phi2tt(comp,i,j,k) = phi2xx + phi2yy + phi2zz - ( 0.5*lambda2*(phi2MagSqr - pow(eta2,2)) + 0.5*lambda12*(phi1MagSqr - pow(eta1,2)) )*phi2(comp,tNow,i,j,k);


                        // Calculate first order derivatives for energy

                        // phix[comp] = ( cos(theta(0,tNow,i,j,k))*phi(comp,tNow,i+1,j,k) + c[comp]*sin(theta(0,tNow,i,j,k))*phi(comp+c[comp],tNow,i+1,j,k) 
                        //              - cos(theta(0,tNow,i-1,j,k))*phi(comp,tNow,i-1,j,k) + c[comp]*sin(theta(0,tNow,i-1,j,k))*phi(comp+c[comp],tNow,i-1,j,k) )/(2*dx);
                        // phiy[comp] = ( cos(theta(1,tNow,i,j,k))*phi(comp,tNow,i,j+1,k) + c[comp]*sin(theta(1,tNow,i,j,k))*phi(comp+c[comp],tNow,i,j+1,k) 
                        //              - cos(theta(1,tNow,i,j-1,k))*phi(comp,tNow,i,j-1,k) + c[comp]*sin(theta(1,tNow,i,j-1,k))*phi(comp+c[comp],tNow,i,j-1,k) )/(2*dy);

                        //Try without central differencing

                        phi1x[comp] = ( phi1(comp,tNow,i,j,k) - cos(q1*theta(0,tNow,i-1,j,k))*phi1(comp,tNow,i-1,j,k) + c[comp]*sin(q1*theta(0,tNow,i-1,j,k))*phi1(comp+c[comp],tNow,i-1,j,k) )/dx;
                        phi1y[comp] = ( phi1(comp,tNow,i,j,k) - cos(q1*theta(1,tNow,i,j-1,k))*phi1(comp,tNow,i,j-1,k) + c[comp]*sin(q1*theta(1,tNow,i,j-1,k))*phi1(comp+c[comp],tNow,i,j-1,k) )/dy;

                        phi2x[comp] = ( phi2(comp,tNow,i,j,k) - cos(q2*theta(0,tNow,i-1,j,k))*phi2(comp,tNow,i-1,j,k) + c[comp]*sin(q2*theta(0,tNow,i-1,j,k))*phi2(comp+c[comp],tNow,i-1,j,k) )/dx;
                        phi2y[comp] = ( phi2(comp,tNow,i,j,k) - cos(q2*theta(1,tNow,i,j-1,k))*phi2(comp,tNow,i,j-1,k) + c[comp]*sin(q2*theta(1,tNow,i,j-1,k))*phi2(comp+c[comp],tNow,i,j-1,k) )/dy;

                        if(k==0){

                            phi1z[comp] = ( phi1(comp,tNow,i,j,k) - cos(q1*theta(2,tNow,i,j,nz-1))*phi1(comp,tNow,i,j,nz-1) + c[comp]*sin(q1*theta(2,tNow,i,j,nz-1))*phi1(comp+c[comp],tNow,i,j,nz-1) )/dz;
                            phi2z[comp] = ( phi2(comp,tNow,i,j,k) - cos(q2*theta(2,tNow,i,j,nz-1))*phi2(comp,tNow,i,j,nz-1) + c[comp]*sin(q2*theta(2,tNow,i,j,nz-1))*phi2(comp+c[comp],tNow,i,j,nz-1) )/dz;

                        } else if(k==nz-1){

                            phi1z[comp] = ( phi1(comp,tNow,i,j,0) - cos(q1*theta(2,tNow,i,j,k-1))*phi1(comp,tNow,i,j,k-1) + c[comp]*sin(q1*theta(2,tNow,i,j,k-1))*phi1(comp+c[comp],tNow,i,j,k-1) )/dz;
                            phi2z[comp] = ( phi2(comp,tNow,i,j,0) - cos(q2*theta(2,tNow,i,j,k-1))*phi2(comp,tNow,i,j,k-1) + c[comp]*sin(q2*theta(2,tNow,i,j,k-1))*phi2(comp+c[comp],tNow,i,j,k-1) )/dz;

                        } else{

                            // phiz[comp] = ( cos(theta(2,tNow,i,j,k))*phi(comp,tNow,i,j,k+1) + c[comp]*sin(theta(2,tNow,i,j,k))*phi(comp+c[comp],tNow,i,j,k+1) 
                            //              - cos(theta(2,tNow,i,j,k-1))*phi(comp,tNow,i,j,k-1) + c[comp]*sin(theta(2,tNow,i,j,k-1))*phi(comp+c[comp],tNow,i,j,k-1) )/(2*dz);

                            phi1z[comp] = ( phi1(comp,tNow,i,j,k) - cos(q1*theta(2,tNow,i,j,k-1))*phi1(comp,tNow,i,j,k-1) + c[comp]*sin(q1*theta(2,tNow,i,j,k-1))*phi1(comp+c[comp],tNow,i,j,k-1) )/dz;
                            phi2z[comp] = ( phi2(comp,tNow,i,j,k) - cos(q2*theta(2,tNow,i,j,k-1))*phi2(comp,tNow,i,j,k-1) + c[comp]*sin(q2*theta(2,tNow,i,j,k-1))*phi2(comp+c[comp],tNow,i,j,k-1) )/dz;

                        }

                        // Updated field value - previous field value to get time derivative with central differencing

                        //phit[comp] = ( 2*phi(comp,tNow,i,j,k) - 2*phi(comp,tPast,i,j,k) + dt*dt*phitt(comp,i,j,k) )/(2*dt);
                        phi1t[comp] = ( phi1(comp,tNow,i,j,k) - phi1(comp,tPast,i,j,k) )/dt;
                        phi2t[comp] = ( phi2(comp,tNow,i,j,k) - phi2(comp,tPast,i,j,k) )/dt;

                    }

                    // Calculate the currents

                    cur1x = q1*( cos(q1*theta(0,tNow,i,j,k))*(phi1(1,tNow,i,j,k)*phi1(0,tNow,i+1,j,k) - phi1(0,tNow,i,j,k)*phi1(1,tNow,i+1,j,k)) +
                                 sin(q1*theta(0,tNow,i,j,k))*(phi1(0,tNow,i,j,k)*phi1(0,tNow,i+1,j,k) + phi1(1,tNow,i,j,k)*phi1(1,tNow,i+1,j,k)) );

                    cur1y = q1*( cos(q1*theta(1,tNow,i,j,k))*(phi1(1,tNow,i,j,k)*phi1(0,tNow,i,j+1,k) - phi1(0,tNow,i,j,k)*phi1(1,tNow,i,j+1,k)) +
                                 sin(q1*theta(1,tNow,i,j,k))*(phi1(0,tNow,i,j,k)*phi1(0,tNow,i,j+1,k) + phi1(1,tNow,i,j,k)*phi1(1,tNow,i,j+1,k)) );

                    cur1z = q1*( cos(q1*theta(2,tNow,i,j,k))*(phi1(1,tNow,i,j,k)*phi1(0,tNow,i,j,k+1) - phi1(0,tNow,i,j,k)*phi1(1,tNow,i,j,k+1)) +
                                 sin(q1*theta(2,tNow,i,j,k))*(phi1(0,tNow,i,j,k)*phi1(0,tNow,i,j,k+1) + phi1(1,tNow,i,j,k)*phi1(1,tNow,i,j,k+1)) );


                    cur2x = q2*( cos(q2*theta(0,tNow,i,j,k))*(phi2(1,tNow,i,j,k)*phi2(0,tNow,i+1,j,k) - phi2(0,tNow,i,j,k)*phi2(1,tNow,i+1,j,k)) +
                                 sin(q2*theta(0,tNow,i,j,k))*(phi2(0,tNow,i,j,k)*phi2(0,tNow,i+1,j,k) + phi2(1,tNow,i,j,k)*phi2(1,tNow,i+1,j,k)) );

                    cur2y = q2*( cos(q2*theta(1,tNow,i,j,k))*(phi2(1,tNow,i,j,k)*phi2(0,tNow,i,j+1,k) - phi2(0,tNow,i,j,k)*phi2(1,tNow,i,j+1,k)) +
                                 sin(q2*theta(1,tNow,i,j,k))*(phi2(0,tNow,i,j,k)*phi2(0,tNow,i,j+1,k) + phi2(1,tNow,i,j,k)*phi2(1,tNow,i,j+1,k)) );

                    cur2z = q2*( cos(q2*theta(2,tNow,i,j,k))*(phi2(1,tNow,i,j,k)*phi2(0,tNow,i,j,k+1) - phi2(0,tNow,i,j,k)*phi2(1,tNow,i,j,k+1)) +
                                 sin(q2*theta(2,tNow,i,j,k))*(phi2(0,tNow,i,j,k)*phi2(0,tNow,i,j,k+1) + phi2(1,tNow,i,j,k)*phi2(1,tNow,i,j,k+1)) );


                    // Calculate the derivatives of the field tensor (lattice version)

                    Fxy_y = ( sin(theta(0,tNow,i,j,k) + theta(1,tNow,i+1,j,k) - theta(0,tNow,i,j+1,k) - theta(1,tNow,i,j,k)) - 
                              sin(theta(0,tNow,i,j-1,k) + theta(1,tNow,i+1,j-1,k) - theta(0,tNow,i,j,k) - theta(1,tNow,i,j-1,k)) )/(dy*dy);

                    Fxz_z = ( sin(theta(0,tNow,i,j,k) + theta(2,tNow,i+1,j,k) - theta(0,tNow,i,j,k+1) - theta(2,tNow,i,j,k)) -
                              sin(theta(0,tNow,i,j,k-1) + theta(2,tNow,i+1,j,k-1) - theta(0,tNow,i,j,k) - theta(2,tNow,i,j,k-1)) )/(dz*dz);

                    Fyx_x = ( sin(theta(1,tNow,i,j,k) + theta(0,tNow,i,j+1,k) - theta(1,tNow,i+1,j,k) - theta(0,tNow,i,j,k)) -
                              sin(theta(1,tNow,i-1,j,k) + theta(0,tNow,i-1,j+1,k) - theta(1,tNow,i,j,k) - theta(0,tNow,i-1,j,k)) )/(dx*dx);

                    Fyz_z = ( sin(theta(1,tNow,i,j,k) + theta(2,tNow,i,j+1,k) - theta(1,tNow,i,j,k+1) - theta(2,tNow,i,j,k)) -
                              sin(theta(1,tNow,i,j,k-1) + theta(2,tNow,i,j+1,k-1) - theta(1,tNow,i,j,k) - theta(2,tNow,i,j,k-1)) )/(dz*dz);

                    Fzx_x = ( sin(theta(2,tNow,i,j,k) + theta(0,tNow,i,j,k+1) - theta(2,tNow,i+1,j,k) - theta(0,tNow,i,j,k)) -
                              sin(theta(2,tNow,i-1,j,k) + theta(0,tNow,i-1,j,k+1) - theta(2,tNow,i,j,k) - theta(0,tNow,i-1,j,k)) )/(dx*dx);

                    Fzy_y = ( sin(theta(2,tNow,i,j,k) + theta(1,tNow,i,j,k+1) - theta(2,tNow,i,j+1,k) - theta(1,tNow,i,j,k)) -
                              sin(theta(2,tNow,i,j-1,k) + theta(1,tNow,i,j-1,k+1) - theta(2,tNow,i,j,k) - theta(1,tNow,i,j-1,k)) )/(dy*dy);

                    thetatt(0,i,j,k) = -2*g*g*(cur1x + cur2x) - Fxy_y - Fxz_z;
                    thetatt(1,i,j,k) = -2*g*g*(cur1y + cur2y) - Fyx_x - Fyz_z;
                    thetatt(2,i,j,k) = -2*g*g*(cur1z + cur2z) - Fzx_x - Fzy_y;

                    // Gauge condition and energy contribution calculations

                    if(g==0){ 

                        gaussDeviation(i,j,k) = 0;
                        thetaDotCont = 0;
                        FCont = 0; 

                    }
                    else{

                        // Calculate deriation from the gauge condition. divTheta/g is the divergence of the gauge field

                        divTheta[tNow] = ( theta(0,tNow,i,j,k) - theta(0,tNow,i-1,j,k) )/(dx*dx) + ( theta(1,tNow,i,j,k) - theta(1,tNow,i,j-1,k) )/(dy*dy) 
                                       + ( theta(2,tNow,i,j,k) - theta(2,tNow,i,j,k-1) )/(dz*dz);
                        divTheta[tPast] = ( theta(0,tPast,i,j,k) - theta(0,tPast,i-1,j,k) )/(dx*dx) + ( theta(1,tPast,i,j,k) - theta(1,tPast,i,j-1,k) )/(dy*dy) 
                                        + ( theta(2,tPast,i,j,k) - theta(2,tPast,i,j,k-1) )/(dz*dz);
                        divThetat = ( divTheta[tNow] - divTheta[tPast] )/dt;

                        gaussDeviation(i,j,k) = divThetat/g - 2*g*( q1*(phi1(1,tNow,i,j,k)*phi1t[0] - phi1(0,tNow,i,j,k)*phi1t[1]) + q2*(phi2(1,tNow,i,j,k)*phi2t[0] - phi2(0,tNow,i,j,k)*phi2t[1]) ); 

                        // Energy contribution terms

                        // These terms are really the time derivatives of the gauge fields (hence factor of grid spacing and g)
                        // thetat[0] = ( 2*theta(0,tNow,i,j,k) - 2*theta(0,tPast,i,j,k) + dt*dt*thetatt(0,i,j,k) )/(2*dt*dx*g);
                        // thetat[1] = ( 2*theta(1,tNow,i,j,k) - 2*theta(1,tPast,i,j,k) + dt*dt*thetatt(1,i,j,k) )/(2*dt*dy*g);
                        // thetat[2] = ( 2*theta(2,tNow,i,j,k) - 2*theta(2,tPast,i,j,k) + dt*dt*thetatt(2,i,j,k) )/(2*dt*dz*g);

                        thetat[0] = ( theta(0,tNow,i,j,k) - theta(0,tPast,i,j,k) )/(g*dt*dx);
                        thetat[1] = ( theta(1,tNow,i,j,k) - theta(1,tPast,i,j,k) )/(g*dt*dy);
                        thetat[2] = ( theta(2,tNow,i,j,k) - theta(2,tPast,i,j,k) )/(g*dt*dz);

                        thetaDotCont = 0.5*( pow(thetat[0],2) + pow(thetat[1],2) + pow(thetat[2],2) ); 

                        // Field strength terms calculated from Wilson loop (factor of half ignored as Fxy=Fyx)

                        Fxy = ( 1 - cos(theta(0,tNow,i,j,k) + theta(1,tNow,i+1,j,k) - theta(0,tNow,i,j+1,k) - theta(1,tNow,i,j,k)) )/pow(g*dx*dy,2);
                        Fxz = ( 1 - cos(theta(0,tNow,i,j,k) + theta(2,tNow,i+1,j,k) - theta(0,tNow,i,j,k+1) - theta(2,tNow,i,j,k)) )/pow(g*dx*dz,2);
                        Fyz = ( 1 - cos(theta(1,tNow,i,j,k) + theta(2,tNow,i,j+1,k) - theta(1,tNow,i,j,k+1) - theta(2,tNow,i,j,k)) )/pow(g*dy*dz,2);

                        FCont = Fxy + Fxz + Fyz;

                    }

                    deviationParameter += abs(gaussDeviation(i,j,k))/((nx-2)*(ny-2)*(nz-2));  // Divide by number of grid points being summed over to get average magnitude of gauss deviation.


                    energydensity(i,j,k) = pow(phi1t[0],2) + pow(phi1t[1],2) + pow(phi1x[0],2) + pow(phi1x[1],2) + pow(phi1y[0],2) + pow(phi1y[1],2) + pow(phi1z[0],2) + pow(phi1z[1],2)
                                         + pow(phi2t[0],2) + pow(phi2t[1],2) + pow(phi2x[0],2) + pow(phi2x[1],2) + pow(phi2y[0],2) + pow(phi2y[1],2) + pow(phi2z[0],2) + pow(phi2z[1],2)
                                         + thetaDotCont + FCont + 0.25*lambda1*pow(phi1MagSqr-pow(eta1,2),2) + 0.25*lambda2*pow(phi2MagSqr-pow(eta2,2),2)
                                         + 0.5*lambda12*(phi1MagSqr-pow(eta1,2))*(phi2MagSqr-pow(eta2,2));

                    energy += dx*dy*dz*energydensity(i,j,k);

                }
            }
        }


        //////////////////////////////////////////////////////////////////////////////////////////
        //                            Absorbing boundary conditions
        //////////////////////////////////////////////////////////////////////////////////////////

        // if(xyBC=="absorbing"){

        //     // Set along x boundary

        //     for(j=1;j<ny-1;j++){
        //         for(k=1;k<nz-1;k++){
        //             for(comp=0;comp<2;comp++){

        //                 // x<0 boundary


        //                 phiyy = ( cos(theta(1,tNow,0,j,k))*phi(comp,tNow,0,j+1,k) + c[comp]*sin(theta(1,tNow,0,j,k))*phi(comp+c[comp],tNow,0,j+1,k) - 2*phi(comp,tNow,0,j,k)
        //                         + cos(theta(1,tNow,0,j-1,k))*phi(comp,tNow,0,j-1,k) - c[comp]*sin(theta(1,tNow,0,j-1,k))*phi(comp+c[comp],tNow,0,j-1,k) )/(dy*dy);

        //                 phizz = ( cos(theta(2,tNow,0,j,k))*phi(comp,tNow,0,j,k+1) + c[comp]*sin(theta(2,tNow,0,j,k))*phi(comp+c[comp],tNow,0,j,k+1) - 2*phi(comp,tNow,0,j,k)
        //                         + cos(theta(2,tNow,0,j,k-1))*phi(comp,tNow,0,j,k-1) - c[comp]*sin(theta(2,tNow,0,j,k-1))*phi(comp+c[comp],tNow,0,j,k-1) )/(dz*dz);

        //                 phitx = ( cos(theta(0,tNow,0,j,k))*phi(comp,tNow,1,j,k) + c[comp]*sin(theta(0,tNow,0,j,k))*phi(comp+c[comp],tNow,1,j,k) - phi(comp,tNow,0,j,k) 
        //                         - cos(theta(0,tPast,0,j,k))*phi(comp,tPast,1,j,k) - c[comp]*sin(theta(0,tPast,0,j,k))*phi(comp+c[comp],tPast,1,j,k) + phi(comp,tPast,0,j,k) )/(dt*dx);

        //                 phitt(comp,0,j,k) = phitx + 0.5*(phiyy + phizz);

        //                 // x>0 boundary

        //                 phiyy = ( cos(theta(1,tNow,nx-1,j,k))*phi(comp,tNow,nx-1,j+1,k) + c[comp]*sin(theta(1,tNow,nx-1,j,k))*phi(comp+c[comp],tNow,nx-1,j+1,k) - 2*phi(comp,tNow,nx-1,j,k) 
        //                         + cos(theta(1,tNow,nx-1,j-1,k))*phi(comp,tNow,nx-1,j-1,k) - c[comp]*sin(theta(1,tNow,nx-1,j-1,k))*phi(comp+c[comp],tNow,nx-1,j-1,k) )/(dy*dy);

        //                 phizz = ( cos(theta(2,tNow,nx-1,j,k))*phi(comp,tNow,nx-1,j,k+1) + c[comp]*sin(theta(2,tNow,nx-1,j,k))*phi(comp+c[comp],tNow,nx-1,j,k+1) - 2*phi(comp,tNow,nx-1,j,k) 
        //                         + cos(theta(2,tNow,nx-1,j,k-1))*phi(comp,tNow,nx-1,j,k-1) - c[comp]*sin(theta(2,tNow,nx-1,j,k-1))*phi(comp+c[comp],tNow,nx-1,j,k-1) )/(dz*dz);

        //                 phitx = ( phi(comp,tNow,nx-1,j,k) - cos(theta(0,tNow,nx-2,j,k))*phi(comp,tNow,nx-2,j,k) + c[comp]*sin(theta(0,tNow,nx-2,j,k))*phi(comp+c[comp],tNow,nx-2,j,k) 
        //                         - phi(comp,tPast,nx-1,j,k) + cos(theta(0,tPast,nx-2,j,k))*phi(comp,tPast,nx-2,j,k) - c[comp]*sin(theta(0,tPast,nx-2,j,k))*phi(comp+c[comp],tPast,nx-2,j,k) )/(dt*dx);

        //                 phitt(comp,nx-1,j,k) = -phitx + 0.5*(phiyy + phizz);

        //             }

        //             for(comp=0;comp<3;comp++){

        //                 // x<0 boundary

        //                 thetayy = ( theta(comp,tNow,0,j+1,k) - 2*theta(comp,tNow,0,j,k) + theta(comp,tNow,0,j-1,k) )/(dy*dy);
        //                 thetazz = ( theta(comp,tNow,0,j,k+1) - 2*theta(comp,tNow,0,j,k) + theta(comp,tNow,0,j,k-1) )/(dz*dz);

        //                 thetatx = ( theta(comp,tNow,1,j,k) - theta(comp,tNow,0,j,k) - theta(comp,tPast,1,j,k) + theta(comp,tPast,0,j,k) )/(dt*dx);

        //                 thetatt(comp,0,j,k) = thetatx + 0.5*(thetayy + thetazz);

        //                 // x>0 boundary

        //                 thetayy = ( theta(comp,tNow,nx-1,j+1,k) - 2*theta(comp,tNow,nx-1,j,k) + theta(comp,tNow,nx-1,j-1,k) )/(dy*dy);
        //                 thetazz = ( theta(comp,tNow,nx-1,j,k+1) - 2*theta(comp,tNow,nx-1,j,k) + theta(comp,tNow,nx-1,j,k-1) )/(dz*dz);

        //                 thetatx = ( theta(comp,tNow,nx-1,j,k) - theta(comp,tNow,nx-2,j,k) - theta(comp,tPast,nx-1,j,k) + theta(comp,tPast,nx-2,j,k) )/(dt*dx);

        //                 thetatt(comp,nx-1,j,k) = -thetatx + 0.5*(thetayy + thetazz);

        //             }

        //         }
        //     }

        //     // Set along y boundary

        //     for(i=1;i<nx-1;i++){
        //         for(k=1;k<nz-1;k++){
        //             for(comp=0;comp<2;comp++){

        //                 // y<0 boundary

        //                 phixx = ( cos(theta(0,tNow,i,0,k))*phi(comp,tNow,i+1,0,k) + c[comp]*sin(theta(0,tNow,i,0,k))*phi(comp+c[comp],tNow,i+1,0,k) - 2*phi(comp,tNow,i,0,k) 
        //                         + cos(theta(0,tNow,i-1,0,k))*phi(comp,tNow,i-1,0,k) - c[comp]*sin(theta(0,tNow,i-1,0,k))*phi(comp+c[comp],tNow,i-1,0,k) )/(dx*dx);

        //                 phizz = ( cos(theta(2,tNow,i,0,k))*phi(comp,tNow,i,0,k+1) + c[comp]*sin(theta(2,tNow,i,0,k))*phi(comp+c[comp],tNow,i,0,k+1) - 2*phi(comp,tNow,i,0,k) 
        //                         + cos(theta(2,tNow,i,0,k-1))*phi(comp,tNow,i,0,k-1) - c[comp]*sin(theta(2,tNow,i,0,k-1))*phi(comp+c[comp],tNow,i,0,k-1) )/(dz*dz);

        //                 phity = ( cos(theta(1,tNow,i,0,k))*phi(comp,tNow,i,1,k) + c[comp]*sin(theta(1,tNow,i,0,k))*phi(comp+c[comp],tNow,i,1,k) - phi(comp,tNow,i,0,k) 
        //                         - cos(theta(1,tPast,i,0,k))*phi(comp,tPast,i,1,k) - c[comp]*sin(theta(1,tPast,i,0,k))*phi(comp+c[comp],tPast,i,1,k) + phi(comp,tPast,i,0,k) )/(dt*dy);

        //                 phitt(comp,i,0,k) = phity + 0.5*(phixx + phizz);


        //                 // y>0 boundary

        //                 phixx = ( cos(theta(0,tNow,i,ny-1,k))*phi(comp,tNow,i+1,ny-1,k) + c[comp]*sin(theta(0,tNow,i,ny-1,k))*phi(comp+c[comp],tNow,i+1,ny-1,k) - 2*phi(comp,tNow,i,ny-1,k) 
        //                         + cos(theta(0,tNow,i-1,ny-1,k))*phi(comp,tNow,i-1,ny-1,k) - c[comp]*sin(theta(0,tNow,i-1,ny-1,k))*phi(comp+c[comp],tNow,i-1,ny-1,k) )/(dx*dx);

        //                 phizz = ( cos(theta(2,tNow,i,ny-1,k))*phi(comp,tNow,i,ny-1,k+1) + c[comp]*sin(theta(2,tNow,i,ny-1,k))*phi(comp+c[comp],tNow,i,ny-1,k+1) - 2*phi(comp,tNow,i,ny-1,k) 
        //                         + cos(theta(2,tNow,i,ny-1,k-1))*phi(comp,tNow,i,ny-1,k-1) - c[comp]*sin(theta(2,tNow,i,ny-1,k-1))*phi(comp+c[comp],tNow,i,ny-1,k-1) )/(dz*dz);

        //                 phity = ( phi(comp,tNow,i,ny-1,k) - cos(theta(1,tNow,i,ny-2,k))*phi(comp,tNow,i,ny-2,k) + c[comp]*sin(theta(1,tNow,i,ny-2,k))*phi(comp+c[comp],tNow,i,ny-2,k) 
        //                         - phi(comp,tPast,i,ny-1,k) + cos(theta(1,tPast,i,ny-2,k))*phi(comp,tPast,i,ny-2,k) - c[comp]*sin(theta(1,tPast,i,ny-2,k))*phi(comp+c[comp],tPast,i,ny-2,k) )/(dt*dy);

        //                 phitt(comp,i,ny-1,k) = -phity + 0.5*(phixx + phizz);

        //             }

        //             for(comp=0;comp<3;comp++){

        //                 // y<0 boundary

        //                 thetaxx = ( theta(comp,tNow,i+1,0,k) - 2*theta(comp,tNow,i,0,k) + theta(comp,tNow,i-1,0,k) )/(dx*dx);
        //                 thetazz = ( theta(comp,tNow,i,0,k+1) - 2*theta(comp,tNow,i,0,k) + theta(comp,tNow,i,0,k-1) )/(dz*dz);

        //                 thetaty = ( theta(comp,tNow,i,1,k) - theta(comp,tNow,i,0,k) - theta(comp,tPast,i,1,k) + theta(comp,tPast,i,0,k) )/(dt*dy);

        //                 thetatt(comp,i,0,k) = thetaty + 0.5*(thetaxx + thetazz);

        //                 // y>0 boundary

        //                 thetaxx = ( theta(comp,tNow,i+1,ny-1,k) - 2*theta(comp,tNow,i,ny-1,k) + theta(comp,tNow,i-1,ny-1,k) )/(dx*dx);
        //                 thetazz = ( theta(comp,tNow,i,ny-1,k+1) - 2*theta(comp,tNow,i,ny-1,k) + theta(comp,tNow,i,ny-1,k-1) )/(dz*dz);

        //                 thetaty = ( theta(comp,tNow,i,ny-1,k) - theta(comp,tNow,i,ny-2,k) - theta(comp,tPast,i,ny-1,k) + theta(comp,tPast,i,ny-2,k) )/(dt*dy);

        //                 thetatt(comp,i,ny-1,k) = -thetaty + 0.5*(thetaxx + thetazz);

        //             }

        //         }
        //     }

        //     // Assign corners. Calculated by adding both absorbing both condition equations and subtracting 1/2 times the wave equation

        //     for(k=1;k<nz-1;k++){
        //         for(comp=0;comp<2;comp++){

        //             // x,y<0 corner

        //             phizz = ( cos(theta(2,tNow,0,0,k))*phi(comp,tNow,0,0,k+1) + c[comp]*sin(theta(2,tNow,0,0,k))*phi(comp+c[comp],tNow,0,0,k+1) - 2*phi(comp,tNow,0,0,k) 
        //                     + cos(theta(2,tNow,0,0,k-1))*phi(comp,tNow,0,0,k-1) - c[comp]*sin(theta(2,tNow,0,0,k-1))*phi(comp+c[comp],tNow,0,0,k-1) )/(dz*dz);

        //             phitx = ( cos(theta(0,tNow,0,0,k))*phi(comp,tNow,1,0,k) + c[comp]*sin(theta(0,tNow,0,0,k))*phi(comp+c[comp],tNow,1,0,k) - phi(comp,tNow,0,0,k) 
        //                     - cos(theta(0,tPast,0,0,k))*phi(comp,tPast,1,0,k) - c[comp]*sin(theta(0,tPast,0,0,k))*phi(comp+c[comp],tPast,1,0,k) + phi(comp,tPast,0,0,k) )/(dt*dx);

        //             phity = ( cos(theta(1,tNow,0,0,k))*phi(comp,tNow,0,1,k) + c[comp]*sin(theta(1,tNow,0,0,k))*phi(comp+c[comp],tNow,0,1,k) - phi(comp,tNow,0,0,k) 
        //                     - cos(theta(1,tPast,0,0,k))*phi(comp,tPast,0,1,k) - c[comp]*sin(theta(1,tPast,0,0,k))*phi(comp+c[comp],tPast,0,1,k) + phi(comp,tPast,0,0,k) )/(dt*dy);

        //             phitt(comp,0,0,k) = ( phizz + 2*(phitx + phity) )/3;

        //             // x<0,y>0 corner

        //             phizz = ( cos(theta(2,tNow,0,ny-1,k))*phi(comp,tNow,0,ny-1,k+1) + c[comp]*sin(theta(2,tNow,0,ny-1,k))*phi(comp+c[comp],tNow,0,ny-1,k+1) - 2*phi(comp,tNow,0,ny-1,k) 
        //                     + cos(theta(2,tNow,0,ny-1,k-1))*phi(comp,tNow,0,ny-1,k-1) - c[comp]*sin(theta(2,tNow,0,ny-1,k-1))*phi(comp+c[comp],tNow,0,ny-1,k-1) )/(dz*dz);

        //             phitx = ( cos(theta(0,tNow,0,ny-1,k))*phi(comp,tNow,1,ny-1,k) + c[comp]*sin(theta(0,tNow,0,ny-1,k))*phi(comp+c[comp],tNow,1,ny-1,k) - phi(comp,tNow,0,ny-1,k) 
        //                     - cos(theta(0,tPast,0,ny-1,k))*phi(comp,tPast,1,ny-1,k) - c[comp]*sin(theta(0,tPast,0,ny-1,k))*phi(comp+c[comp],tPast,1,ny-1,k) + phi(comp,tPast,0,ny-1,k) )/(dt*dx);

        //             phity = ( phi(comp,tNow,0,ny-1,k) - cos(theta(1,tNow,0,ny-2,k))*phi(comp,tNow,0,ny-2,k) + c[comp]*sin(theta(1,tNow,0,ny-2,k))*phi(comp+c[comp],tNow,0,ny-2,k) 
        //                     - phi(comp,tPast,0,ny-1,k) + cos(theta(1,tPast,0,ny-2,k))*phi(comp,tPast,0,ny-2,k) - c[comp]*sin(theta(1,tPast,0,ny-2,k))*phi(comp+c[comp],tPast,0,ny-2,k) )/(dt*dy);

        //             phitt(comp,0,ny-1,k) = ( phizz + 2*(phitx - phity) )/3;

        //             // x>0,y<0 corner

        //             phizz = ( cos(theta(2,tNow,nx-1,0,k))*phi(comp,tNow,nx-1,0,k+1) + c[comp]*sin(theta(2,tNow,nx-1,0,k))*phi(comp+c[comp],tNow,nx-1,0,k+1) - 2*phi(comp,tNow,nx-1,0,k) 
        //                     + cos(theta(2,tNow,nx-1,0,k-1))*phi(comp,tNow,nx-1,0,k-1) - c[comp]*sin(theta(2,tNow,nx-1,0,k-1))*phi(comp+c[comp],tNow,nx-1,0,k-1) )/(dz*dz);

        //             phitx = ( phi(comp,tNow,nx-1,0,k) - cos(theta(0,tNow,nx-2,0,k))*phi(comp,tNow,nx-2,0,k) + c[comp]*sin(theta(0,tNow,nx-2,0,k))*phi(comp+c[comp],tNow,nx-2,0,k) 
        //                     - phi(comp,tPast,nx-1,0,k) + cos(theta(0,tPast,nx-2,0,k))*phi(comp,tPast,nx-2,0,k) - c[comp]*sin(theta(0,tPast,nx-2,0,k))*phi(comp+c[comp],tPast,nx-2,0,k) )/(dt*dx);

        //             phity = ( cos(theta(1,tNow,nx-1,0,k))*phi(comp,tNow,nx-1,1,k) + c[comp]*sin(theta(1,tNow,nx-1,0,k))*phi(comp+c[comp],tNow,nx-1,1,k) - phi(comp,tNow,nx-1,0,k) 
        //                     - cos(theta(1,tPast,nx-1,0,k))*phi(comp,tPast,nx-1,1,k) - c[comp]*sin(theta(1,tPast,nx-1,0,k))*phi(comp+c[comp],tPast,nx-1,1,k) + phi(comp,tPast,nx-1,0,k) )/(dt*dy);

        //             phitt(comp,nx-1,0,k) = ( phizz - 2*(phitx - phity) )/3;

        //             // x,y>0 corner

        //             phizz = ( cos(theta(2,tNow,nx-1,ny-1,k))*phi(comp,tNow,nx-1,ny-1,k+1) + c[comp]*sin(theta(2,tNow,nx-1,ny-1,k))*phi(comp+c[comp],tNow,nx-1,ny-1,k+1) - 2*phi(comp,tNow,nx-1,ny-1,k) 
        //                     + cos(theta(2,tNow,nx-1,ny-1,k-1))*phi(comp,tNow,nx-1,ny-1,k-1) - c[comp]*sin(theta(2,tNow,nx-1,ny-1,k-1))*phi(comp+c[comp],tNow,nx-1,ny-1,k-1) )/(dz*dz);

        //             phitx = ( phi(comp,tNow,nx-1,ny-1,k) - cos(theta(0,tNow,nx-2,ny-1,k))*phi(comp,tNow,nx-2,ny-1,k) + c[comp]*sin(theta(0,tNow,nx-2,ny-1,k))*phi(comp+c[comp],tNow,nx-2,ny-1,k) 
        //                     - phi(comp,tPast,nx-1,ny-1,k) + cos(theta(0,tPast,nx-1,ny-1,k))*phi(comp,tPast,nx-2,ny-1,k) - c[comp]*sin(theta(0,tPast,nx-2,ny-1,k))*phi(comp+c[comp],tPast,nx-2,ny-1,k)
        //                     )/(dt*dx);

        //             phity = ( phi(comp,tNow,nx-1,ny-1,k) - cos(theta(1,tNow,nx-1,ny-2,k))*phi(comp,tNow,nx-1,ny-2,k) + c[comp]*sin(theta(1,tNow,nx-1,ny-2,k))*phi(comp+c[comp],tNow,nx-1,ny-2,k) 
        //                     - phi(comp,tPast,nx-1,ny-1,k) + cos(theta(1,tPast,nx-1,ny-2,k))*phi(comp,tPast,nx-1,ny-2,k) - c[comp]*sin(theta(1,tPast,nx-1,ny-2,k))*phi(comp+c[comp],tPast,nx-1,ny-2,k) 
        //                     )/(dt*dy);

        //             phitt(comp,nx-1,ny-1,k) = ( phizz - 2*(phitx + phity) )/3;

        //         }

        //         for(comp=0;comp<3;comp++){

        //             // x,y<0 corner

        //             thetazz = ( theta(comp,tNow,0,0,k+1) - 2*theta(comp,tNow,0,0,k) + theta(comp,tNow,0,0,k-1) )/(dz*dz);

        //             thetatx = ( theta(comp,tNow,1,0,k) - theta(comp,tNow,0,0,k) - theta(comp,tPast,1,0,k) + theta(comp,tPast,0,0,k) )/(dt*dx);
        //             thetaty = ( theta(comp,tNow,0,1,k) - theta(comp,tNow,0,0,k) - theta(comp,tPast,0,1,k) + theta(comp,tPast,0,0,k) )/(dt*dy);

        //             thetatt(comp,0,0,k) = ( thetazz + 2*(thetatx + thetaty) )/3;

        //             // x<0,y>0 corner

        //             thetazz = ( theta(comp,tNow,0,ny-1,k+1) - 2*theta(comp,tNow,0,ny-1,k) + theta(comp,tNow,0,ny-1,k-1) )/(dz*dz);

        //             thetatx = ( theta(comp,tNow,1,ny-1,k) - theta(comp,tNow,0,ny-1,k) - theta(comp,tPast,1,ny-1,k) + theta(comp,tPast,0,ny-1,k) )/(dt*dx);
        //             thetaty = ( theta(comp,tNow,0,ny-1,k) - theta(comp,tNow,0,ny-2,k) - theta(comp,tPast,0,ny-1,k) + theta(comp,tPast,0,ny-2,k) )/(dt*dy);

        //             thetatt(comp,0,ny-1,k) = ( thetazz + 2*(thetatx - thetaty) )/3;

        //             // x>0,y<0 corner

        //             thetazz = ( theta(comp,tNow,nx-1,0,k+1) - 2*theta(comp,tNow,nx-1,0,k) + theta(comp,tNow,nx-1,0,k-1) )/(dz*dz);

        //             thetatx = ( theta(comp,tNow,nx-1,0,k) - theta(comp,tNow,nx-2,0,k) - theta(comp,tPast,nx-1,0,k) + theta(comp,tPast,nx-2,0,k) )/(dt*dx);
        //             thetaty = ( theta(comp,tNow,nx-1,1,k) - theta(comp,tNow,nx-1,0,k) - theta(comp,tPast,nx-1,1,k) + theta(comp,tPast,nx-1,0,k) )/(dt*dy);

        //             thetatt(comp,nx-1,0,k) = ( thetazz - 2*(thetatx - thetaty) )/3;

        //             // x,y>0 corner

        //             thetazz = ( theta(comp,tNow,nx-1,ny-1,k+1) - 2*theta(comp,tNow,nx-1,ny-1,k) + theta(comp,tNow,nx-1,ny-1,k-1) )/(dz*dz);

        //             thetatx = ( theta(comp,tNow,nx-1,ny-1,k) - theta(comp,tNow,nx-2,ny-1,k) - theta(comp,tPast,nx-1,ny-1,k) + theta(comp,tPast,nx-2,ny-1,k) )/(dt*dx);
        //             thetaty = ( theta(comp,tNow,nx-1,ny-1,k) - theta(comp,tNow,nx-1,ny-2,k) - theta(comp,tPast,nx-1,ny-1,k) + theta(comp,tPast,nx-1,ny-2,k) )/(dt*dy);

        //             thetatt(comp,nx-1,ny-1,k) = ( thetazz - 2*(thetatx + thetaty) )/3;

        //         }

        //     }


        // }



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

            #pragma omp parallel for default(none) shared(phi1,phi2,theta,phi1tt,phi2tt,thetatt,tNow,tPast) private(j,k,comp)
            for(i=0;i<nx;i++){
                for(j=0;j<ny;j++){
                    for(k=1;k<nz-1;k++){

                        for(comp=0;comp<2;comp++){

                            phi1(comp,tPast,i,j,k) = 2*phi1(comp,tNow,i,j,k) - phi1(comp,tPast,i,j,k) + dt*dt*phi1tt(comp,i,j,k);
                            phi2(comp,tPast,i,j,k) = 2*phi2(comp,tNow,i,j,k) - phi2(comp,tPast,i,j,k) + dt*dt*phi2tt(comp,i,j,k);

                        }

                        for(comp=0;comp<3;comp++){

                            theta(comp,tPast,i,j,k) = 2*theta(comp,tNow,i,j,k) - theta(comp,tPast,i,j,k) + dt*dt*thetatt(comp,i,j,k);

                        }
                    }
                }
            }

        } else{

            #pragma omp parallel for default(none) shared(phi1,phi2,theta,phi1tt,phi2tt,thetatt,tNow,tPast) private(j,k,comp)
            for(i=1;i<nx-1;i++){
                for(j=1;j<ny-1;j++){
                    for(k=1;k<nz-1;k++){

                        for(comp=0;comp<2;comp++){

                            phi1(comp,tPast,i,j,k) = 2*phi1(comp,tNow,i,j,k) - phi1(comp,tPast,i,j,k) + dt*dt*phi1tt(comp,i,j,k);
                            phi2(comp,tPast,i,j,k) = 2*phi2(comp,tNow,i,j,k) - phi2(comp,tPast,i,j,k) + dt*dt*phi2tt(comp,i,j,k);

                        }

                        for(comp=0;comp<3;comp++){

                            theta(comp,tPast,i,j,k) = 2*theta(comp,tNow,i,j,k) - theta(comp,tPast,i,j,k) + dt*dt*thetatt(comp,i,j,k);

                        }

                    }
                }
            }

        }


        valsPerLoop << energy << " " << deviationParameter/((nx-2)*(ny-2)*(nz-2)) << endl;


        if(makeGif && TimeStep%saveFreq == 0){

            ss.str(string());
            ss << gifFrame;
            string gifDataPath = dir_path + "/GifData/gifData_" + ss.str() + ".txt";
            ofstream gifData (gifDataPath.c_str());
            gifFrame+=1;

            for(i=0;i<nx;i++){
                for(j=0;j<ny;j++){
                    for(k=0;k<nz;k++){

                        gifData << pow(phi1(0,tPast,i,j,k),2) + pow(phi1(1,tPast,i,j,k),2) << " " << pow(phi2(0,tPast,i,j,k),2) + pow(phi2(1,tPast,i,j,k),2) << endl;

                    }
                }
            }

        }

    }

    cout << "\rTimestep " << nt << " completed." << endl;

    for(i=0;i<nx;i++){
        for(j=0;j<ny;j++){
            for(k=0;k<nz;k++){

                finalField << phi1(0,tPast,i,j,k) << " " << phi1(1,tPast,i,j,k) << " " << phi2(0,tPast,i,j,k) << " " << phi2(1,tPast,i,j,k) << " " 
                           << theta(0,tPast,i,j,k) << " " << theta(1,tPast,i,j,k) << " " << theta(2,tPast,i,j,k) << endl;
                test1 << energydensity(i,j,k) << endl;

            }
        }
    }


    gettimeofday(&end,NULL);

    cout << "Time taken: " << end.tv_sec - start.tv_sec << "s" << endl;


	return 0;

}


