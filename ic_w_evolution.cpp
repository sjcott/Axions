#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <sys/time.h>
#include <omp.h>
#include <complex>
#include <vector>

#include "array.hpp"
#include "ic_func.hpp"

using namespace std;

// Perturb a straight, global string solution, evolve the system and investigate how it radiates.

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                 		  Parameters & Declarations
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Some parts of code may assume nx,ny and nz are odd numbers

const int nx = 201;
const int ny = 201;
const int nz = 201;
const int nt = 2001;
const double dx = 0.35;
const double dy = 0.35;
const double dz = 0.35;
const double dt = 0.05;

const double lambda = 1;
const double eta = 1;
const double g = 0;

const int damped_nt = 0; // Number of time steps for which damping is imposed. Useful for random initial conditions
const double dampFac = 0; // magnitude of damping term, unclear how strong to make this

const bool makeGif = true;
const int saveFreq = 10;
const int countRate = 10;

// Centaurus doesn't like me putting these in openmp shared() if they are const and my laptop/telesto don't like it if I dont.
// My solution has been to just remove the const declaration
string xyBC = "periodic"; // Allows for "neumann" (covariant derivatives set to zero), "absorbing", "periodic" or fixed (any other string will choose this option) boundary conditions.
string zBC = "periodic"; // Allows for "neumann" (covariant derivatives set to zero), "periodic" or fixed (any other string will choose this option) boundary conditions.

const bool stringPos = false;
const bool stringDetect = true; // true if you want the code to find which faces the strings pass through. May try to calculate string length later too.
const bool splitLength = true;


int main(){

	struct timeval start, end;
    gettimeofday(&start, NULL);

    string file_path = __FILE__;
    string dir_path = file_path.substr(0,file_path.find_last_of('/'));
    stringstream ss;

    string input;
    //cout << "Enter a tag for output files: " << flush;
    //cin >> input;
    input = "test";

    string finalFieldPath = dir_path + "/Data/finalField.txt";
    string valsPerLoopPath = dir_path + "/Data/valsPerLoop_" + input + ".txt";
    //string valsPerLoopPath = dir_path + "/Data/valsPerLoop_test.txt";
    string test1Path = dir_path + "/test1.txt";
    string test2Path = dir_path + "/test2.txt";
    string powerdensityOutPath = dir_path + "/Data/powerdensityOut.txt";

    ofstream finalField (finalFieldPath.c_str());
    ofstream valsPerLoop (valsPerLoopPath.c_str());
    ofstream test1 (test1Path.c_str());
    ofstream test2 (test2Path.c_str());
    ofstream powerdensityOut (powerdensityOutPath.c_str());

    twoArray ic = ic_func();
    Array phi(2,2,nx,ny,nz,0.0), theta(3,2,nx,ny,nz,0.0);
    phi = ic.array1;
    theta = ic.array2;

    // // Use a scope block to deallocate the memory used by the twoArray struct "ic" once data has been transferred to the Arrays 
    // {

    //     twoArray ic = ic_func(); // Declare struct and receive the two Arrays phi and theta from the initial conditions function
    //     phi = ic.array1;
    //     theta = ic.array2;

    // }

    Array phitt(2,nx,ny,nz,0.0), thetatt(3,nx,ny,nz,0.0), energydensity(2,nx,ny,nz,0.0), powerdensity(nx,ny,nz,0.0), gaussDeviation(nx,ny,nz,0.0);
    int comp, i, j, k, TimeStep, gifFrame, tNow, tPast, s, counter, x0, y0, z0, xEdge1, xEdge2, yEdge1, yEdge2, zEdge1, zEdge2, im, jm, km, ip, jp, kp;
    double phixx, phiyy, phizz, phiMagSqr, phix[2], phiy[2], phiz[2], curx, cury, curz, Fxy_y, Fxz_z, Fyx_x, Fyz_z, Fzx_x, Fzy_y, phit[2], energy, phitx, phity, thetat[3], divTheta[2], divThetat,
           thetaDotCont, Fxy, Fxz, Fyz, FCont, deviationParameter, thetaxx, thetayy, thetazz, thetatx, thetaty, damp, stringLength[2];
    int c[2] = {1,-1}; // Useful definition to allow covariant deviative to be calculated when looping over components.

    x0 = int(0.5*(nx-1));
    y0 = int(0.5*(ny-1));
    z0 = int(0.5*(nz-1));

    gifFrame = 0;
    counter = 0;

    for(TimeStep=0;TimeStep<nt;TimeStep++){

        //if((1000*TimeStep)%(nt-1)==0){
        if(TimeStep>counter){

            cout << "\rTimestep " << TimeStep-1 << " completed." << flush;

            counter += countRate;

        }

        // Is damping switched on or not?
        if(TimeStep<damped_nt){ damp = dampFac; }
        else{ damp = 0; }

        tNow = (TimeStep+1)%2;
        tPast = TimeStep%2;

        if(xyBC == "periodic"){ 

            xEdge1 = 0;
            xEdge2 = nx;
            yEdge1 = 0;
            yEdge2 = ny; 

        } else{

            xEdge1 = 1;
            xEdge2 = nx-1;
            yEdge1 = 1;
            yEdge2 = ny-1;

        }

        if(zBC == "periodic"){

            zEdge1 = 0;
            zEdge2 = nz;

        } else{

            zEdge1 = 1;
            zEdge2 = nz-1;

        }

        // Set boundary conditions. Different loop sizes so that corner allocation is done correctly.

        // Sets Neumann boundary conditions on x and y

        if(xyBC=="neumann"){

            for(i=1;i<nx-1;i++){
                for(k=zEdge1;k<zEdge2;k++){

                    // Doing theta first because it is used in one of the covariant derivative of phi
                    for(comp=0;comp<3;comp++){

                        theta(comp,tNow,i,0,k) = theta(comp,tNow,i,1,k);
                        theta(comp,tNow,i,ny-1,k) = theta(comp,tNow,i,ny-2,k);

                    }

                    for(comp=0;comp<2;comp++){

                        phi(comp,tNow,i,0,k) = cos(theta(1,tNow,i,0,k))*phi(comp,tNow,i,1,k) + c[comp]*sin(theta(1,tNow,i,0,k))*phi(comp+c[comp],tNow,i,1,k);
                        phi(comp,tNow,i,ny-1,k) = cos(theta(1,tNow,i,ny-2,k))*phi(comp,tNow,i,ny-2,k) - c[comp]*sin(theta(1,tNow,i,ny-2,k))*phi(comp+c[comp],tNow,i,ny-2,k);

                    }

                }
            }

            for(j=0;j<ny;j++){
                for(k=zEdge1;k<zEdge2;k++){

                    // Doing theta first because it is used in one of the covariant derivative of phi
                    for(comp=0;comp<3;comp++){

                        theta(comp,tNow,0,j,k) = theta(comp,tNow,1,j,k);
                        theta(comp,tNow,nx-1,j,k) = theta(comp,tNow,nx-2,j,k);

                    }

                    for(comp=0;comp<2;comp++){

                        phi(comp,tNow,0,j,k) = cos(theta(0,tNow,0,j,k))*phi(comp,tNow,1,j,k) + c[comp]*sin(theta(0,tNow,0,j,k))*phi(comp+c[comp],tNow,1,j,k);
                        phi(comp,tNow,nx-1,j,k) = cos(theta(0,tNow,nx-2,j,k))*phi(comp,tNow,nx-2,j,k) - c[comp]*sin(theta(0,tNow,nx-2,j,k))*phi(comp+c[comp],tNow,nx-2,j,k);

                    }

                }
            }

        }

        // Set Neumann derivatives at the z boundaries

        if(zBC == "neumann"){

            for(i=0;i<nz;i++){
                for(j=0;j<ny;j++){

                    // Do theta first as it's used in covariant derivatives of phi
                    for(comp=0;comp<3;comp++){

                        theta(comp,tNow,i,j,0) = theta(comp,tNow,i,j,1);
                        theta(comp,tNow,i,j,nz-1) = theta(comp,tNow,i,j,nz-2);

                    }

                    for(comp=0;comp<2;comp++){

                        phi(comp,tNow,i,j,0) = cos(theta(2,tNow,i,j,0))*phi(comp,tNow,i,j,1) + c[comp]*sin(theta(2,tNow,i,j,0))*phi(comp+c[comp],tNow,i,j,1);
                        phi(comp,tNow,i,j,nz-1) = cos(theta(2,tNow,i,j,nz-2))*phi(comp,tNow,i,j,nz-2) - c[comp]*sin(theta(2,tNow,i,j,nz-2))*phi(comp+c[comp],tNow,i,j,nz-2);

                    }

                }
            }

        }


        // Calculate time derivatives using EoMs

        energy = 0;
        deviationParameter = 0;

        #pragma omp parallel for reduction(+:energy,deviationParameter) default(none) shared(phi,theta,phitt,thetatt,energydensity,powerdensity,gaussDeviation,tNow,tPast,c,TimeStep,damp,xEdge1, \
        xEdge2,yEdge1,yEdge2,zEdge1,zEdge2) \
        private(phiMagSqr,phixx,phiyy,phizz,j,k,comp,phix,phiy,phiz,phit,curx,cury,curz,Fxy_y,Fxz_z,Fyx_x,Fyz_z,Fzx_x,Fzy_y,divTheta,divThetat,thetat,thetaDotCont,Fxy,Fxz,Fyz,FCont,im,jm,km,ip,jp,kp)

        for(i=xEdge1;i<xEdge2;i++){

            // im and ip are used to index the neigbouring points (and same for jm,jp,km,km)

            ip = (i+1)%nx;
            im = (i-1+nx)%nx;

            for(j=yEdge1;j<yEdge2;j++){

                jp = (j+1)%ny;
                jm = (j-1+ny)%ny;

                for(k=zEdge1;k<zEdge2;k++){

                    kp = (k+1)%nz;
                    km = (k-1+nz)%nz;

                    phiMagSqr = pow(phi(0,tNow,i,j,k),2) + pow(phi(1,tNow,i,j,k),2);

                    // Loop over phi components

                    for(comp=0;comp<2;comp++){

                         // 2nd order spatial derivatives calculated with 2nd order finite difference

                        // c is 1 when comp = 0 and -1 when comp = 1

                        phixx = ( cos(theta(0,tNow,i,j,k))*phi(comp,tNow,ip,j,k) + c[comp]*sin(theta(0,tNow,i,j,k))*phi(comp+c[comp],tNow,ip,j,k) - 2*phi(comp,tNow,i,j,k)
                                + cos(theta(0,tNow,im,j,k))*phi(comp,tNow,im,j,k) - c[comp]*sin(theta(0,tNow,im,j,k))*phi(comp+c[comp],tNow,im,j,k) )/(dx*dx);

                        phiyy = ( cos(theta(1,tNow,i,j,k))*phi(comp,tNow,i,jp,k) + c[comp]*sin(theta(1,tNow,i,j,k))*phi(comp+c[comp],tNow,i,jp,k) - 2*phi(comp,tNow,i,j,k)
                                + cos(theta(1,tNow,i,jm,k))*phi(comp,tNow,i,jm,k) - c[comp]*sin(theta(1,tNow,i,jm,k))*phi(comp+c[comp],tNow,i,jm,k) )/(dy*dy);
                        
                        phizz = ( cos(theta(2,tNow,i,j,k))*phi(comp,tNow,i,j,kp) + c[comp]*sin(theta(2,tNow,i,j,k))*phi(comp+c[comp],tNow,i,j,kp) - 2*phi(comp,tNow,i,j,k) 
                                + cos(theta(2,tNow,i,j,km))*phi(comp,tNow,i,j,km) - c[comp]*sin(theta(2,tNow,i,j,km))*phi(comp+c[comp],tNow,i,j,km) )/(dz*dz);


                        phit[comp] = ( phi(comp,tNow,i,j,k) - phi(comp,tPast,i,j,k) )/dt;

                        // Calculate second order time derivatives

                        phitt(comp,i,j,k) = phixx + phiyy + phizz - 0.5*lambda*(phiMagSqr - pow(eta,2))*phi(comp,tNow,i,j,k) - damp*phit[comp];


                        // Calculate first order derivatives for energy

                        // phix[comp] = ( cos(theta(0,tNow,i,j,k))*phi(comp,tNow,i+1,j,k) + c[comp]*sin(theta(0,tNow,i,j,k))*phi(comp+c[comp],tNow,i+1,j,k) 
                        //              - cos(theta(0,tNow,i-1,j,k))*phi(comp,tNow,i-1,j,k) + c[comp]*sin(theta(0,tNow,i-1,j,k))*phi(comp+c[comp],tNow,i-1,j,k) )/(2*dx);
                        // phiy[comp] = ( cos(theta(1,tNow,i,j,k))*phi(comp,tNow,i,j+1,k) + c[comp]*sin(theta(1,tNow,i,j,k))*phi(comp+c[comp],tNow,i,j+1,k) 
                        //              - cos(theta(1,tNow,i,j-1,k))*phi(comp,tNow,i,j-1,k) + c[comp]*sin(theta(1,tNow,i,j-1,k))*phi(comp+c[comp],tNow,i,j-1,k) )/(2*dy);

                        //Try without central differencing

                        phix[comp] = ( phi(comp,tNow,i,j,k) - cos(theta(0,tNow,im,j,k))*phi(comp,tNow,im,j,k) + c[comp]*sin(theta(0,tNow,im,j,k))*phi(comp+c[comp],tNow,im,j,k) )/dx;
                        phiy[comp] = ( phi(comp,tNow,i,j,k) - cos(theta(1,tNow,i,jm,k))*phi(comp,tNow,i,jm,k) + c[comp]*sin(theta(1,tNow,i,jm,k))*phi(comp+c[comp],tNow,i,jm,k) )/dy;
                        phiz[comp] = ( phi(comp,tNow,i,j,k) - cos(theta(2,tNow,i,j,km))*phi(comp,tNow,i,j,km) + c[comp]*sin(theta(2,tNow,i,j,km))*phi(comp+c[comp],tNow,i,j,km) )/dz;


                    }

                    // Calculate the currents

                    curx = cos(theta(0,tNow,i,j,k))*(phi(1,tNow,i,j,k)*phi(0,tNow,ip,j,k) - phi(0,tNow,i,j,k)*phi(1,tNow,ip,j,k)) +
                           sin(theta(0,tNow,i,j,k))*(phi(0,tNow,i,j,k)*phi(0,tNow,ip,j,k) + phi(1,tNow,i,j,k)*phi(1,tNow,ip,j,k));

                    cury = cos(theta(1,tNow,i,j,k))*(phi(1,tNow,i,j,k)*phi(0,tNow,i,jp,k) - phi(0,tNow,i,j,k)*phi(1,tNow,i,jp,k)) +
                           sin(theta(1,tNow,i,j,k))*(phi(0,tNow,i,j,k)*phi(0,tNow,i,jp,k) + phi(1,tNow,i,j,k)*phi(1,tNow,i,jp,k));

                    curz = cos(theta(2,tNow,i,j,k))*(phi(1,tNow,i,j,k)*phi(0,tNow,i,j,kp) - phi(0,tNow,i,j,k)*phi(1,tNow,i,j,kp)) +
                           sin(theta(2,tNow,i,j,k))*(phi(0,tNow,i,j,k)*phi(0,tNow,i,j,kp) + phi(1,tNow,i,j,k)*phi(1,tNow,i,j,kp));


                    // Calculate the derivatives of the field tensor (lattice version)

                    Fxy_y = ( sin(theta(0,tNow,i,j,k) + theta(1,tNow,ip,j,k) - theta(0,tNow,i,jp,k) - theta(1,tNow,i,j,k)) - 
                              sin(theta(0,tNow,i,jm,k) + theta(1,tNow,ip,jm,k) - theta(0,tNow,i,j,k) - theta(1,tNow,i,jm,k)) )/(dy*dy);

                    Fxz_z = ( sin(theta(0,tNow,i,j,k) + theta(2,tNow,ip,j,k) - theta(0,tNow,i,j,kp) - theta(2,tNow,i,j,k)) -
                              sin(theta(0,tNow,i,j,km) + theta(2,tNow,ip,j,km) - theta(0,tNow,i,j,k) - theta(2,tNow,i,j,km)) )/(dz*dz);

                    Fyx_x = ( sin(theta(1,tNow,i,j,k) + theta(0,tNow,i,jp,k) - theta(1,tNow,ip,j,k) - theta(0,tNow,i,j,k)) -
                              sin(theta(1,tNow,im,j,k) + theta(0,tNow,im,jp,k) - theta(1,tNow,i,j,k) - theta(0,tNow,im,j,k)) )/(dx*dx);

                    Fyz_z = ( sin(theta(1,tNow,i,j,k) + theta(2,tNow,i,jp,k) - theta(1,tNow,i,j,kp) - theta(2,tNow,i,j,k)) -
                              sin(theta(1,tNow,i,j,km) + theta(2,tNow,i,jp,km) - theta(1,tNow,i,j,k) - theta(2,tNow,i,j,km)) )/(dz*dz);

                    Fzx_x = ( sin(theta(2,tNow,i,j,k) + theta(0,tNow,i,j,kp) - theta(2,tNow,ip,j,k) - theta(0,tNow,i,j,k)) -
                              sin(theta(2,tNow,im,j,k) + theta(0,tNow,im,j,kp) - theta(2,tNow,i,j,k) - theta(0,tNow,im,j,k)) )/(dx*dx);

                    Fzy_y = ( sin(theta(2,tNow,i,j,k) + theta(1,tNow,i,j,kp) - theta(2,tNow,i,jp,k) - theta(1,tNow,i,j,k)) -
                              sin(theta(2,tNow,i,jm,k) + theta(1,tNow,i,jm,kp) - theta(2,tNow,i,j,k) - theta(1,tNow,i,jm,k)) )/(dy*dy);

                    thetatt(0,i,j,k) = -2*g*g*curx - Fxy_y - Fxz_z - damp*( theta(0,tNow,i,j,k) - theta(0,tPast,i,j,k) )/dt;
                    thetatt(1,i,j,k) = -2*g*g*cury - Fyx_x - Fyz_z - damp*( theta(1,tNow,i,j,k) - theta(1,tPast,i,j,k) )/dt;
                    thetatt(2,i,j,k) = -2*g*g*curz - Fzx_x - Fzy_y - damp*( theta(2,tNow,i,j,k) - theta(2,tPast,i,j,k) )/dt;

                    // Gauge condition and energy contribution calculations

                    if(g==0){ 

                        gaussDeviation(i,j,k) = 0;
                        thetaDotCont = 0;
                        FCont = 0; 

                    }
                    else{

                        // Calculate deriation from the gauge condition. divTheta/g is the divergence of the gauge field

                        divTheta[tNow] = ( theta(0,tNow,i,j,k) - theta(0,tNow,im,j,k) )/(dx*dx) + ( theta(1,tNow,i,j,k) - theta(1,tNow,i,jm,k) )/(dy*dy) 
                                       + ( theta(2,tNow,i,j,k) - theta(2,tNow,i,j,km) )/(dz*dz);
                        divTheta[tPast] = ( theta(0,tPast,i,j,k) - theta(0,tPast,im,j,k) )/(dx*dx) + ( theta(1,tPast,i,j,k) - theta(1,tPast,i,jm,k) )/(dy*dy) 
                                        + ( theta(2,tPast,i,j,k) - theta(2,tPast,i,j,km) )/(dz*dz);
                        divThetat = ( divTheta[tNow] - divTheta[tPast] )/dt;

                        gaussDeviation(i,j,k) = divThetat/g - 2*g*(phi(1,tNow,i,j,k)*phit[0] - phi(0,tNow,i,j,k)*phit[1]); 

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

                        Fxy = ( 1 - cos(theta(0,tNow,i,j,k) + theta(1,tNow,ip,j,k) - theta(0,tNow,i,jp,k) - theta(1,tNow,i,j,k)) )/pow(g*dx*dy,2);
                        Fxz = ( 1 - cos(theta(0,tNow,i,j,k) + theta(2,tNow,ip,j,k) - theta(0,tNow,i,j,kp) - theta(2,tNow,i,j,k)) )/pow(g*dx*dz,2);
                        Fyz = ( 1 - cos(theta(1,tNow,i,j,k) + theta(2,tNow,i,jp,k) - theta(1,tNow,i,j,kp) - theta(2,tNow,i,j,k)) )/pow(g*dy*dz,2);

                        FCont = Fxy + Fxz + Fyz;

                    }

                    deviationParameter += abs(gaussDeviation(i,j,k));


                    energydensity(tNow,i,j,k) = pow(phit[0],2) + pow(phit[1],2) + pow(phix[0],2) + pow(phix[1],2) + pow(phiy[0],2) + pow(phiy[1],2) + pow(phiz[0],2) + pow(phiz[1],2)
                                              + thetaDotCont + FCont + 0.25*lambda*pow(phiMagSqr-pow(eta,2),2);

                    powerdensity(i,j,k) = (energydensity(tNow,i,j,k) - energydensity(tPast,i,j,k))/dt;

                    energy += dx*dy*dz*energydensity(tNow,i,j,k);

                }
            }
        }

        //////////////////////////////////////////////////////////////////////////////////////////
        //                            Absorbing boundary conditions
        //////////////////////////////////////////////////////////////////////////////////////////

        if(xyBC=="absorbing"){

            // Set along x boundary

            for(j=1;j<ny-1;j++){
                for(k=zEdge1;k<zEdge2;k++){

                    kp = (k+1)%nz;
                    km = (k-1+nz)%nz;

                    for(comp=0;comp<2;comp++){



                        // x<0 boundary


                        phiyy = ( cos(theta(1,tNow,0,j,k))*phi(comp,tNow,0,j+1,k) + c[comp]*sin(theta(1,tNow,0,j,k))*phi(comp+c[comp],tNow,0,j+1,k) - 2*phi(comp,tNow,0,j,k)
                                + cos(theta(1,tNow,0,j-1,k))*phi(comp,tNow,0,j-1,k) - c[comp]*sin(theta(1,tNow,0,j-1,k))*phi(comp+c[comp],tNow,0,j-1,k) )/(dy*dy);

                        phizz = ( cos(theta(2,tNow,0,j,k))*phi(comp,tNow,0,j,kp) + c[comp]*sin(theta(2,tNow,0,j,k))*phi(comp+c[comp],tNow,0,j,kp) - 2*phi(comp,tNow,0,j,k)
                                + cos(theta(2,tNow,0,j,km))*phi(comp,tNow,0,j,km) - c[comp]*sin(theta(2,tNow,0,j,km))*phi(comp+c[comp],tNow,0,j,km) )/(dz*dz);

                        phitx = ( cos(theta(0,tNow,0,j,k))*phi(comp,tNow,1,j,k) + c[comp]*sin(theta(0,tNow,0,j,k))*phi(comp+c[comp],tNow,1,j,k) - phi(comp,tNow,0,j,k) 
                                - cos(theta(0,tPast,0,j,k))*phi(comp,tPast,1,j,k) - c[comp]*sin(theta(0,tPast,0,j,k))*phi(comp+c[comp],tPast,1,j,k) + phi(comp,tPast,0,j,k) )/(dt*dx);

                        phitt(comp,0,j,k) = phitx + 0.5*(phiyy + phizz);

                        // x>0 boundary

                        phiyy = ( cos(theta(1,tNow,nx-1,j,k))*phi(comp,tNow,nx-1,j+1,k) + c[comp]*sin(theta(1,tNow,nx-1,j,k))*phi(comp+c[comp],tNow,nx-1,j+1,k) - 2*phi(comp,tNow,nx-1,j,k) 
                                + cos(theta(1,tNow,nx-1,j-1,k))*phi(comp,tNow,nx-1,j-1,k) - c[comp]*sin(theta(1,tNow,nx-1,j-1,k))*phi(comp+c[comp],tNow,nx-1,j-1,k) )/(dy*dy);

                        phizz = ( cos(theta(2,tNow,nx-1,j,k))*phi(comp,tNow,nx-1,j,kp) + c[comp]*sin(theta(2,tNow,nx-1,j,k))*phi(comp+c[comp],tNow,nx-1,j,kp) - 2*phi(comp,tNow,nx-1,j,k) 
                                + cos(theta(2,tNow,nx-1,j,km))*phi(comp,tNow,nx-1,j,km) - c[comp]*sin(theta(2,tNow,nx-1,j,km))*phi(comp+c[comp],tNow,nx-1,j,km) )/(dz*dz);

                        phitx = ( phi(comp,tNow,nx-1,j,k) - cos(theta(0,tNow,nx-2,j,k))*phi(comp,tNow,nx-2,j,k) + c[comp]*sin(theta(0,tNow,nx-2,j,k))*phi(comp+c[comp],tNow,nx-2,j,k) 
                                - phi(comp,tPast,nx-1,j,k) + cos(theta(0,tPast,nx-2,j,k))*phi(comp,tPast,nx-2,j,k) - c[comp]*sin(theta(0,tPast,nx-2,j,k))*phi(comp+c[comp],tPast,nx-2,j,k) )/(dt*dx);

                        phitt(comp,nx-1,j,k) = -phitx + 0.5*(phiyy + phizz);

                    }

                    for(comp=0;comp<3;comp++){

                        // x<0 boundary

                        thetayy = ( theta(comp,tNow,0,j+1,k) - 2*theta(comp,tNow,0,j,k) + theta(comp,tNow,0,j-1,k) )/(dy*dy);
                        thetazz = ( theta(comp,tNow,0,j,kp) - 2*theta(comp,tNow,0,j,k) + theta(comp,tNow,0,j,km) )/(dz*dz);

                        thetatx = ( theta(comp,tNow,1,j,k) - theta(comp,tNow,0,j,k) - theta(comp,tPast,1,j,k) + theta(comp,tPast,0,j,k) )/(dt*dx);

                        thetatt(comp,0,j,k) = thetatx + 0.5*(thetayy + thetazz);

                        // x>0 boundary

                        thetayy = ( theta(comp,tNow,nx-1,j+1,k) - 2*theta(comp,tNow,nx-1,j,k) + theta(comp,tNow,nx-1,j-1,k) )/(dy*dy);
                        thetazz = ( theta(comp,tNow,nx-1,j,kp) - 2*theta(comp,tNow,nx-1,j,k) + theta(comp,tNow,nx-1,j,km) )/(dz*dz);

                        thetatx = ( theta(comp,tNow,nx-1,j,k) - theta(comp,tNow,nx-2,j,k) - theta(comp,tPast,nx-1,j,k) + theta(comp,tPast,nx-2,j,k) )/(dt*dx);

                        thetatt(comp,nx-1,j,k) = -thetatx + 0.5*(thetayy + thetazz);

                    }

                }
            }

            // Set along y boundary

            for(i=1;i<nx-1;i++){
                for(k=zEdge1;k<zEdge2;k++){

                    kp = (k+1)%nz;
                    km = (k-1+nz)%nz;

                    for(comp=0;comp<2;comp++){

                        // y<0 boundary

                        phixx = ( cos(theta(0,tNow,i,0,k))*phi(comp,tNow,i+1,0,k) + c[comp]*sin(theta(0,tNow,i,0,k))*phi(comp+c[comp],tNow,i+1,0,k) - 2*phi(comp,tNow,i,0,k) 
                                + cos(theta(0,tNow,i-1,0,k))*phi(comp,tNow,i-1,0,k) - c[comp]*sin(theta(0,tNow,i-1,0,k))*phi(comp+c[comp],tNow,i-1,0,k) )/(dx*dx);

                        phizz = ( cos(theta(2,tNow,i,0,k))*phi(comp,tNow,i,0,kp) + c[comp]*sin(theta(2,tNow,i,0,k))*phi(comp+c[comp],tNow,i,0,kp) - 2*phi(comp,tNow,i,0,k) 
                                + cos(theta(2,tNow,i,0,km))*phi(comp,tNow,i,0,km) - c[comp]*sin(theta(2,tNow,i,0,km))*phi(comp+c[comp],tNow,i,0,km) )/(dz*dz);

                        phity = ( cos(theta(1,tNow,i,0,k))*phi(comp,tNow,i,1,k) + c[comp]*sin(theta(1,tNow,i,0,k))*phi(comp+c[comp],tNow,i,1,k) - phi(comp,tNow,i,0,k) 
                                - cos(theta(1,tPast,i,0,k))*phi(comp,tPast,i,1,k) - c[comp]*sin(theta(1,tPast,i,0,k))*phi(comp+c[comp],tPast,i,1,k) + phi(comp,tPast,i,0,k) )/(dt*dy);

                        phitt(comp,i,0,k) = phity + 0.5*(phixx + phizz);


                        // y>0 boundary

                        phixx = ( cos(theta(0,tNow,i,ny-1,k))*phi(comp,tNow,i+1,ny-1,k) + c[comp]*sin(theta(0,tNow,i,ny-1,k))*phi(comp+c[comp],tNow,i+1,ny-1,k) - 2*phi(comp,tNow,i,ny-1,k) 
                                + cos(theta(0,tNow,i-1,ny-1,k))*phi(comp,tNow,i-1,ny-1,k) - c[comp]*sin(theta(0,tNow,i-1,ny-1,k))*phi(comp+c[comp],tNow,i-1,ny-1,k) )/(dx*dx);

                        phizz = ( cos(theta(2,tNow,i,ny-1,k))*phi(comp,tNow,i,ny-1,kp) + c[comp]*sin(theta(2,tNow,i,ny-1,k))*phi(comp+c[comp],tNow,i,ny-1,kp) - 2*phi(comp,tNow,i,ny-1,k) 
                                + cos(theta(2,tNow,i,ny-1,km))*phi(comp,tNow,i,ny-1,km) - c[comp]*sin(theta(2,tNow,i,ny-1,km))*phi(comp+c[comp],tNow,i,ny-1,km) )/(dz*dz);

                        phity = ( phi(comp,tNow,i,ny-1,k) - cos(theta(1,tNow,i,ny-2,k))*phi(comp,tNow,i,ny-2,k) + c[comp]*sin(theta(1,tNow,i,ny-2,k))*phi(comp+c[comp],tNow,i,ny-2,k) 
                                - phi(comp,tPast,i,ny-1,k) + cos(theta(1,tPast,i,ny-2,k))*phi(comp,tPast,i,ny-2,k) - c[comp]*sin(theta(1,tPast,i,ny-2,k))*phi(comp+c[comp],tPast,i,ny-2,k) )/(dt*dy);

                        phitt(comp,i,ny-1,k) = -phity + 0.5*(phixx + phizz);

                    }

                    for(comp=0;comp<3;comp++){

                        // y<0 boundary

                        thetaxx = ( theta(comp,tNow,i+1,0,k) - 2*theta(comp,tNow,i,0,k) + theta(comp,tNow,i-1,0,k) )/(dx*dx);
                        thetazz = ( theta(comp,tNow,i,0,kp) - 2*theta(comp,tNow,i,0,k) + theta(comp,tNow,i,0,km) )/(dz*dz);

                        thetaty = ( theta(comp,tNow,i,1,k) - theta(comp,tNow,i,0,k) - theta(comp,tPast,i,1,k) + theta(comp,tPast,i,0,k) )/(dt*dy);

                        thetatt(comp,i,0,k) = thetaty + 0.5*(thetaxx + thetazz);

                        // y>0 boundary

                        thetaxx = ( theta(comp,tNow,i+1,ny-1,k) - 2*theta(comp,tNow,i,ny-1,k) + theta(comp,tNow,i-1,ny-1,k) )/(dx*dx);
                        thetazz = ( theta(comp,tNow,i,ny-1,kp) - 2*theta(comp,tNow,i,ny-1,k) + theta(comp,tNow,i,ny-1,km) )/(dz*dz);

                        thetaty = ( theta(comp,tNow,i,ny-1,k) - theta(comp,tNow,i,ny-2,k) - theta(comp,tPast,i,ny-1,k) + theta(comp,tPast,i,ny-2,k) )/(dt*dy);

                        thetatt(comp,i,ny-1,k) = -thetaty + 0.5*(thetaxx + thetazz);

                    }

                }
            }

            // Assign corners. Calculated by adding both absorbing both condition equations and subtracting 1/2 times the wave equation

            for(k=zEdge1;k<zEdge2;k++){

                kp = (k+1)%nz;
                km = (k-1+nz)%nz;

                for(comp=0;comp<2;comp++){

                    // x,y<0 corner

                    phizz = ( cos(theta(2,tNow,0,0,k))*phi(comp,tNow,0,0,kp) + c[comp]*sin(theta(2,tNow,0,0,k))*phi(comp+c[comp],tNow,0,0,kp) - 2*phi(comp,tNow,0,0,k) 
                            + cos(theta(2,tNow,0,0,km))*phi(comp,tNow,0,0,km) - c[comp]*sin(theta(2,tNow,0,0,km))*phi(comp+c[comp],tNow,0,0,km) )/(dz*dz);

                    phitx = ( cos(theta(0,tNow,0,0,k))*phi(comp,tNow,1,0,k) + c[comp]*sin(theta(0,tNow,0,0,k))*phi(comp+c[comp],tNow,1,0,k) - phi(comp,tNow,0,0,k) 
                            - cos(theta(0,tPast,0,0,k))*phi(comp,tPast,1,0,k) - c[comp]*sin(theta(0,tPast,0,0,k))*phi(comp+c[comp],tPast,1,0,k) + phi(comp,tPast,0,0,k) )/(dt*dx);

                    phity = ( cos(theta(1,tNow,0,0,k))*phi(comp,tNow,0,1,k) + c[comp]*sin(theta(1,tNow,0,0,k))*phi(comp+c[comp],tNow,0,1,k) - phi(comp,tNow,0,0,k) 
                            - cos(theta(1,tPast,0,0,k))*phi(comp,tPast,0,1,k) - c[comp]*sin(theta(1,tPast,0,0,k))*phi(comp+c[comp],tPast,0,1,k) + phi(comp,tPast,0,0,k) )/(dt*dy);

                    phitt(comp,0,0,k) = ( phizz + 2*(phitx + phity) )/3;

                    // x<0,y>0 corner

                    phizz = ( cos(theta(2,tNow,0,ny-1,k))*phi(comp,tNow,0,ny-1,kp) + c[comp]*sin(theta(2,tNow,0,ny-1,k))*phi(comp+c[comp],tNow,0,ny-1,kp) - 2*phi(comp,tNow,0,ny-1,k) 
                            + cos(theta(2,tNow,0,ny-1,km))*phi(comp,tNow,0,ny-1,km) - c[comp]*sin(theta(2,tNow,0,ny-1,km))*phi(comp+c[comp],tNow,0,ny-1,km) )/(dz*dz);

                    phitx = ( cos(theta(0,tNow,0,ny-1,k))*phi(comp,tNow,1,ny-1,k) + c[comp]*sin(theta(0,tNow,0,ny-1,k))*phi(comp+c[comp],tNow,1,ny-1,k) - phi(comp,tNow,0,ny-1,k) 
                            - cos(theta(0,tPast,0,ny-1,k))*phi(comp,tPast,1,ny-1,k) - c[comp]*sin(theta(0,tPast,0,ny-1,k))*phi(comp+c[comp],tPast,1,ny-1,k) + phi(comp,tPast,0,ny-1,k) )/(dt*dx);

                    phity = ( phi(comp,tNow,0,ny-1,k) - cos(theta(1,tNow,0,ny-2,k))*phi(comp,tNow,0,ny-2,k) + c[comp]*sin(theta(1,tNow,0,ny-2,k))*phi(comp+c[comp],tNow,0,ny-2,k) 
                            - phi(comp,tPast,0,ny-1,k) + cos(theta(1,tPast,0,ny-2,k))*phi(comp,tPast,0,ny-2,k) - c[comp]*sin(theta(1,tPast,0,ny-2,k))*phi(comp+c[comp],tPast,0,ny-2,k) )/(dt*dy);

                    phitt(comp,0,ny-1,k) = ( phizz + 2*(phitx - phity) )/3;

                    // x>0,y<0 corner

                    phizz = ( cos(theta(2,tNow,nx-1,0,k))*phi(comp,tNow,nx-1,0,kp) + c[comp]*sin(theta(2,tNow,nx-1,0,k))*phi(comp+c[comp],tNow,nx-1,0,kp) - 2*phi(comp,tNow,nx-1,0,k) 
                            + cos(theta(2,tNow,nx-1,0,km))*phi(comp,tNow,nx-1,0,km) - c[comp]*sin(theta(2,tNow,nx-1,0,km))*phi(comp+c[comp],tNow,nx-1,0,km) )/(dz*dz);

                    phitx = ( phi(comp,tNow,nx-1,0,k) - cos(theta(0,tNow,nx-2,0,k))*phi(comp,tNow,nx-2,0,k) + c[comp]*sin(theta(0,tNow,nx-2,0,k))*phi(comp+c[comp],tNow,nx-2,0,k) 
                            - phi(comp,tPast,nx-1,0,k) + cos(theta(0,tPast,nx-2,0,k))*phi(comp,tPast,nx-2,0,k) - c[comp]*sin(theta(0,tPast,nx-2,0,k))*phi(comp+c[comp],tPast,nx-2,0,k) )/(dt*dx);

                    phity = ( cos(theta(1,tNow,nx-1,0,k))*phi(comp,tNow,nx-1,1,k) + c[comp]*sin(theta(1,tNow,nx-1,0,k))*phi(comp+c[comp],tNow,nx-1,1,k) - phi(comp,tNow,nx-1,0,k) 
                            - cos(theta(1,tPast,nx-1,0,k))*phi(comp,tPast,nx-1,1,k) - c[comp]*sin(theta(1,tPast,nx-1,0,k))*phi(comp+c[comp],tPast,nx-1,1,k) + phi(comp,tPast,nx-1,0,k) )/(dt*dy);

                    phitt(comp,nx-1,0,k) = ( phizz - 2*(phitx - phity) )/3;

                    // x,y>0 corner

                    phizz = ( cos(theta(2,tNow,nx-1,ny-1,k))*phi(comp,tNow,nx-1,ny-1,kp) + c[comp]*sin(theta(2,tNow,nx-1,ny-1,k))*phi(comp+c[comp],tNow,nx-1,ny-1,kp) - 2*phi(comp,tNow,nx-1,ny-1,k) 
                            + cos(theta(2,tNow,nx-1,ny-1,km))*phi(comp,tNow,nx-1,ny-1,km) - c[comp]*sin(theta(2,tNow,nx-1,ny-1,km))*phi(comp+c[comp],tNow,nx-1,ny-1,km) )/(dz*dz);

                    phitx = ( phi(comp,tNow,nx-1,ny-1,k) - cos(theta(0,tNow,nx-2,ny-1,k))*phi(comp,tNow,nx-2,ny-1,k) + c[comp]*sin(theta(0,tNow,nx-2,ny-1,k))*phi(comp+c[comp],tNow,nx-2,ny-1,k) 
                            - phi(comp,tPast,nx-1,ny-1,k) + cos(theta(0,tPast,nx-1,ny-1,k))*phi(comp,tPast,nx-2,ny-1,k) - c[comp]*sin(theta(0,tPast,nx-2,ny-1,k))*phi(comp+c[comp],tPast,nx-2,ny-1,k)
                            )/(dt*dx);

                    phity = ( phi(comp,tNow,nx-1,ny-1,k) - cos(theta(1,tNow,nx-1,ny-2,k))*phi(comp,tNow,nx-1,ny-2,k) + c[comp]*sin(theta(1,tNow,nx-1,ny-2,k))*phi(comp+c[comp],tNow,nx-1,ny-2,k) 
                            - phi(comp,tPast,nx-1,ny-1,k) + cos(theta(1,tPast,nx-1,ny-2,k))*phi(comp,tPast,nx-1,ny-2,k) - c[comp]*sin(theta(1,tPast,nx-1,ny-2,k))*phi(comp+c[comp],tPast,nx-1,ny-2,k) 
                            )/(dt*dy);

                    phitt(comp,nx-1,ny-1,k) = ( phizz - 2*(phitx + phity) )/3;

                }

                for(comp=0;comp<3;comp++){

                    // x,y<0 corner

                    thetazz = ( theta(comp,tNow,0,0,kp) - 2*theta(comp,tNow,0,0,k) + theta(comp,tNow,0,0,km) )/(dz*dz);

                    thetatx = ( theta(comp,tNow,1,0,k) - theta(comp,tNow,0,0,k) - theta(comp,tPast,1,0,k) + theta(comp,tPast,0,0,k) )/(dt*dx);
                    thetaty = ( theta(comp,tNow,0,1,k) - theta(comp,tNow,0,0,k) - theta(comp,tPast,0,1,k) + theta(comp,tPast,0,0,k) )/(dt*dy);

                    thetatt(comp,0,0,k) = ( thetazz + 2*(thetatx + thetaty) )/3;

                    // x<0,y>0 corner

                    thetazz = ( theta(comp,tNow,0,ny-1,kp) - 2*theta(comp,tNow,0,ny-1,k) + theta(comp,tNow,0,ny-1,km) )/(dz*dz);

                    thetatx = ( theta(comp,tNow,1,ny-1,k) - theta(comp,tNow,0,ny-1,k) - theta(comp,tPast,1,ny-1,k) + theta(comp,tPast,0,ny-1,k) )/(dt*dx);
                    thetaty = ( theta(comp,tNow,0,ny-1,k) - theta(comp,tNow,0,ny-2,k) - theta(comp,tPast,0,ny-1,k) + theta(comp,tPast,0,ny-2,k) )/(dt*dy);

                    thetatt(comp,0,ny-1,k) = ( thetazz + 2*(thetatx - thetaty) )/3;

                    // x>0,y<0 corner

                    thetazz = ( theta(comp,tNow,nx-1,0,kp) - 2*theta(comp,tNow,nx-1,0,k) + theta(comp,tNow,nx-1,0,km) )/(dz*dz);

                    thetatx = ( theta(comp,tNow,nx-1,0,k) - theta(comp,tNow,nx-2,0,k) - theta(comp,tPast,nx-1,0,k) + theta(comp,tPast,nx-2,0,k) )/(dt*dx);
                    thetaty = ( theta(comp,tNow,nx-1,1,k) - theta(comp,tNow,nx-1,0,k) - theta(comp,tPast,nx-1,1,k) + theta(comp,tPast,nx-1,0,k) )/(dt*dy);

                    thetatt(comp,nx-1,0,k) = ( thetazz - 2*(thetatx - thetaty) )/3;

                    // x,y>0 corner

                    thetazz = ( theta(comp,tNow,nx-1,ny-1,kp) - 2*theta(comp,tNow,nx-1,ny-1,k) + theta(comp,tNow,nx-1,ny-1,km) )/(dz*dz);

                    thetatx = ( theta(comp,tNow,nx-1,ny-1,k) - theta(comp,tNow,nx-2,ny-1,k) - theta(comp,tPast,nx-1,ny-1,k) + theta(comp,tPast,nx-2,ny-1,k) )/(dt*dx);
                    thetaty = ( theta(comp,tNow,nx-1,ny-1,k) - theta(comp,tNow,nx-1,ny-2,k) - theta(comp,tPast,nx-1,ny-1,k) + theta(comp,tPast,nx-1,ny-2,k) )/(dt*dy);

                    thetatt(comp,nx-1,ny-1,k) = ( thetazz - 2*(thetatx + thetaty) )/3;

                }

            }


        }

        //////////////////////////////////////// Absorbing boundary conditions end ////////////////////////////////////////////////////////////////////////////////////////

        // Find the position of the string along y=z=0.

        double pos;

        if(stringPos){

            double phiSmallest = 1;
            int iString;
            for(i=0;i<nx;i++){

                phiMagSqr = pow(phi(0,tNow,i,y0,z0),2) + pow(phi(1,tNow,i,y0,z0),2);

                if(phiMagSqr<phiSmallest){

                    phiSmallest = phiMagSqr;
                    iString = i;

                }

            }

            double phi_l = pow(phi(0,tNow,iString-1,y0,z0),2) + pow(phi(1,tNow,iString-1,y0,z0),2);
            double phi_r = pow(phi(0,tNow,iString+1,y0,z0),2) + pow(phi(1,tNow,iString+1,y0,z0),2);

            pos = (iString-x0)*dx + dx*( phi_l - phi_r )/( 2*phi_l - 4*phiSmallest + 2*phi_r );

        } else{ pos = 0; }

        // Find where strings intercept grid faces

        if(stringDetect){

            vector<double> xString, yString, zString; // Declares the vectors and clears them at every loop

            #pragma omp parallel default(none) shared(xString,yString,zString,xEdge1,xEdge2,yEdge1,yEdge2,zEdge1,zEdge2,x0,y0,z0,phi,tNow) private(j,k,ip,jp,kp,comp)
            {

                vector<double> xString_private, yString_private, zString_private; // Declares private vectors that will be collected into the shared vectors in a critical section
                double x, y, z, a, b, c, coeff1[2], coeff2[2], coeff3[2], coeff4[2], sol1, sol2;

                #pragma omp for
                for(i=xEdge1;i<xEdge2;i++){

                    ip = (i+1)%nx;
                    x = (i-x0)*dx;

                    for(j=yEdge1;j<yEdge2;j++){

                        jp = (j+1)%ny;
                        y = (j-y0)*dy;

                        for(k=zEdge1;k<zEdge2;k++){

                            kp = (k+1)%nz;
                            z = (k-z0)*dz;

                            // All faces will be checked on the positive side (i.e for z face, look at the face x to x+dx and y to y+dx)
                            // The other three faces will be checked by other points in the grid

                            // Start with the z directed face

                            // Invert a bilinear interpolation (setting Re(phi) and Im(phi) to zero) to get a pair of equations for the position of the zero. Solve that pair of equations.

                            for(comp=0;comp<2;comp++){

                                // Do the same process for the real and imaginary components

                                coeff1[comp] = phi(comp,tNow,ip,jp,k) - phi(comp,tNow,ip,j,k) - phi(comp,tNow,i,jp,k) + phi(comp,tNow,i,j,k);
                                coeff2[comp] = (y+dy)*(phi(comp,tNow,ip,j,k) - phi(comp,tNow,i,j,k)) + y*(phi(comp,tNow,i,jp,k) - phi(comp,tNow,ip,jp,k));
                                coeff3[comp] = (x+dx)*(phi(comp,tNow,i,jp,k) - phi(comp,tNow,i,j,k)) + x*(phi(comp,tNow,ip,j,k) - phi(comp,tNow,ip,jp,k));
                                coeff4[comp] = (x+dx)*( (y+dy)*phi(comp,tNow,i,j,k) - y*phi(comp,tNow,i,jp,k) ) - x*( (y+dy)*phi(comp,tNow,ip,j,k) - y*phi(comp,tNow,ip,jp,k) );

                            }

                            // Substituting one equation into the other gives a quadratic equation for one of the coordinates. Now calculate the coefficients of the quadratic.

                            a = coeff1[1]*coeff3[0] - coeff1[0]*coeff3[1];  
                            b = coeff1[1]*coeff4[0] + coeff2[1]*coeff3[0] - coeff1[0]*coeff4[1] - coeff2[0]*coeff3[1];
                            c = coeff2[1]*coeff4[0] - coeff2[0]*coeff4[1];

                            // Check if a=0, if so the equation is simpler than quadratic --> sol2 (y in this case) is just -c/b unless b is zero as well
                            if(a==0){

                                if(b!=0){

                                    // There is just one solution

                                    sol2 = -c/b; // y location of string

                                    // Does this lie inside the face?
                                    if(sol2>=y && sol2<=y+dy){

                                        sol1 = -(coeff3[0]*sol2 + coeff4[0])/(coeff1[0]*sol2 + coeff2[0]); // x location of string

                                        // Does this lie inside the face?
                                        if(sol1>=x && sol1<=x+dx){

                                            // The string intersects the face so add it to the vectors

                                            xString_private.push_back(sol1);
                                            yString_private.push_back(sol2);
                                            zString_private.push_back(z);

                                        }

                                    }

                                }

                                // Otherwise no solutions

                            } else if(b*b - 4*a*c >= 0){

                                // There will be two solutions (or a repeated, this may cause trouble). First solution is

                                sol2 = ( -b + sqrt(b*b - 4*a*c) )/(2*a);

                                if(sol2>=y && sol2<=y+dy){

                                    sol1 = -(coeff3[0]*sol2 + coeff4[0])/(coeff1[0]*sol2 + coeff2[0]);

                                    if(sol1>=x && sol1<=x+dx){

                                        xString_private.push_back(sol1);
                                        yString_private.push_back(sol2);
                                        zString_private.push_back(z);

                                    }

                                }

                                // Second solution is

                                sol2 = ( -b - sqrt(b*b - 4*a*c) )/(2*a);

                                if(sol2>=y && sol2<=y+dy){

                                    sol1 = -(coeff3[0]*sol2 + coeff4[0])/(coeff1[0]*sol2 + coeff2[0]);

                                    if(sol1>=x && sol1<=x+dx){

                                        xString_private.push_back(sol1);
                                        yString_private.push_back(sol2);
                                        zString_private.push_back(z);

                                    }

                                }

                            }

                            // Otherwise there will only be complex solutions so the string doesn't intersect the face



                            // Now repeat this process for the y directed face

                            for(comp=0;comp<2;comp++){

                                coeff1[comp] = phi(comp,tNow,ip,j,kp) - phi(comp,tNow,ip,j,k) - phi(comp,tNow,i,j,kp) + phi(comp,tNow,i,j,k);
                                coeff2[comp] = (z+dz)*(phi(comp,tNow,ip,j,k) - phi(comp,tNow,i,j,k)) + z*(phi(comp,tNow,i,j,kp) - phi(comp,tNow,ip,j,kp));
                                coeff3[comp] = (x+dx)*(phi(comp,tNow,i,j,kp) - phi(comp,tNow,i,j,k)) + x*(phi(comp,tNow,ip,j,k) - phi(comp,tNow,ip,j,kp));
                                coeff4[comp] = (x+dx)*( (z+dz)*phi(comp,tNow,i,j,k) - z*phi(comp,tNow,i,j,kp) ) - x*( (z+dz)*phi(comp,tNow,ip,j,k) - z*phi(comp,tNow,ip,j,kp) );

                            }

                            a = coeff1[1]*coeff3[0] - coeff1[0]*coeff3[1];
                            b = coeff1[1]*coeff4[0] + coeff2[1]*coeff3[0] - coeff1[0]*coeff4[1] - coeff2[0]*coeff3[1];
                            c = coeff2[1]*coeff4[0] - coeff2[0]*coeff4[1];

                            if(a==0){

                                if(b!=0){

                                    sol2 = -c/b;

                                    if(sol2>=z && sol2<=z+dz){

                                        sol1 = -(coeff3[0]*sol2 + coeff4[0])/(coeff1[0]*sol2 + coeff2[0]);

                                        if(sol1>=x && sol1<=x+dx){

                                            xString_private.push_back(sol1);
                                            yString_private.push_back(y);
                                            zString_private.push_back(sol2);

                                        }

                                    }

                                }

                            } else if(b*b - 4*a*c >= 0){

                                sol2 = ( -b + sqrt(b*b - 4*a*c) )/(2*a);

                                if(sol2>=z && sol2<=z+dz){

                                    sol1 = -(coeff3[0]*sol2 + coeff4[0])/(coeff1[0]*sol2 + coeff2[0]);

                                    if(sol1>=x && sol1<=x+dx){

                                        xString_private.push_back(sol1);
                                        yString_private.push_back(y);
                                        zString_private.push_back(sol2);

                                    }

                                }

                                sol2 = ( -b - sqrt(b*b - 4*a*c) )/(2*a);

                                if(sol2>=z && sol2<=z+dz){

                                    sol1 = -(coeff3[0]*sol2 + coeff4[0])/(coeff1[0]*sol2 + coeff2[0]);

                                    if(sol1>=x && sol1<=x+dx){

                                        xString_private.push_back(sol1);
                                        yString_private.push_back(y);
                                        zString_private.push_back(sol2);

                                    }

                                }

                            }



                            // Now repeat one more time for the x directed face

                            for(comp=0;comp<2;comp++){

                                coeff1[comp] = phi(comp,tNow,i,jp,kp) - phi(comp,tNow,i,jp,k) - phi(comp,tNow,i,j,kp) + phi(comp,tNow,i,j,k);
                                coeff2[comp] = (z+dz)*(phi(comp,tNow,i,jp,k) - phi(comp,tNow,i,j,k)) + z*(phi(comp,tNow,i,j,kp) - phi(comp,tNow,i,jp,kp));
                                coeff3[comp] = (y+dy)*(phi(comp,tNow,i,j,kp) - phi(comp,tNow,i,j,k)) + y*(phi(comp,tNow,i,jp,k) - phi(comp,tNow,i,jp,kp));
                                coeff4[comp] = (y+dy)*( (z+dz)*phi(comp,tNow,i,j,k) - z*phi(comp,tNow,i,j,kp) ) - y*( (z+dz)*phi(comp,tNow,i,jp,k) - z*phi(comp,tNow,i,jp,kp) );

                            }

                            a = coeff1[1]*coeff3[0] - coeff1[0]*coeff3[1];
                            b = coeff1[1]*coeff4[0] + coeff2[1]*coeff3[0] - coeff1[0]*coeff4[1] - coeff2[0]*coeff3[1];
                            c = coeff2[1]*coeff4[0] - coeff2[0]*coeff4[1];

                            if(a==0){

                                if(b!=0){

                                    sol2 = -c/b;

                                    if(sol2>=z && sol2<=z+dz){

                                        sol1 = -(coeff3[0]*sol2 + coeff4[0])/(coeff1[0]*sol2 + coeff2[0]);

                                        if(sol1>=y && sol1<=y+dy){

                                            xString_private.push_back(x);
                                            yString_private.push_back(sol1);
                                            zString_private.push_back(sol2);

                                        }

                                    }

                                }

                            } else if(b*b - 4*a*c >= 0){

                                sol2 = ( -b + sqrt(b*b - 4*a*c) )/(2*a);

                                if(sol2>=z && sol2<=z+dz){

                                    sol1 = -(coeff3[0]*sol2 + coeff4[0])/(coeff1[0]*sol2 + coeff2[0]);

                                    if(sol1>=y && sol1<=y+dy){

                                        xString_private.push_back(x);
                                        yString_private.push_back(sol1);
                                        zString_private.push_back(sol2);

                                    }

                                }

                                sol2 = ( -b - sqrt(b*b - 4*a*c) )/(2*a);

                                if(sol2>=z && sol2<=z+dz){

                                    sol1 = -(coeff3[0]*sol2 + coeff4[0])/(coeff1[0]*sol2 + coeff2[0]);

                                    if(sol1>=y && sol1<=y+dy){

                                        xString_private.push_back(x);
                                        yString_private.push_back(sol1);
                                        zString_private.push_back(sol2);

                                    }

                                }

                            }

                        }
                    }
                }

                // Critical section to merge the detected intersections of each process into one
                #pragma omp critical
                {

                    xString.insert(xString.end(), xString_private.begin(), xString_private.end());
                    yString.insert(yString.end(), yString_private.begin(), yString_private.end());
                    zString.insert(zString.end(), zString_private.begin(), zString_private.end());

                }

            }

            // Now want to calculate the length of string by connecting these points together
            // Need to account for the periodic boundary conditions.
            // Find nearest neighbour to current intersection point and follow it along. If distance > sqrt(dx^2 + dy^2 + dz^2) then it does not count as a neighbour (to avoid linking disconnected strings)

            stringLength[0] = 0;
            stringLength[1] = 0;
            int fInd = 0; // Index of the intersection point we are currently following
            int initialSize = xString.size();
            double beginCoords[3];
            if(initialSize!=0){ beginCoords[0] = xString[fInd]; beginCoords[1] = yString[fInd]; beginCoords[2] = zString[fInd]; }
            for(i=0;i<initialSize;i++){

                bool NextIntersectFound = false;
                double minDistSqr = pow(dx*dx + dy*dy + dz*dz,2);
                int indClosest;

                #pragma omp parallel default(none) shared(NextIntersectFound,minDistSqr,indClosest,xString,yString,zString,fInd,xyBC,zBC)
                {

                    double minDistSqr_private = pow(dx*dx + dy*dy + dz*dz,2);
                    int indClosest_private;

                    #pragma omp for
                    for(j=0;j<xString.size();j++){
                        if(j!=fInd){

                            double distSqr = 0;

                            if(abs(xString[fInd] - xString[j]) <= 0.5*nx*dx){ distSqr = pow(xString[fInd] - xString[j],2); }
                            else if(xyBC == "periodic"){ distSqr = pow(nx*dx - abs(xString[fInd] - xString[j]),2); } // Accounts for periodicity, if using periodic boundaries

                            if(abs(yString[fInd] - yString[j]) <= 0.5*ny*dy){ distSqr += pow(yString[fInd] - yString[j],2); }
                            else if(xyBC == "periodic"){ distSqr += pow(ny*dy - abs(yString[fInd] - yString[j]),2); }

                            if(abs(zString[fInd] - zString[j]) <= 0.5*nz*dz){ distSqr += pow(zString[fInd] - zString[j],2); }
                            else if(zBC == "periodic"){ distSqr += pow(nz*dz - abs(zString[fInd] - zString[j]),2); }

                            // Now check if this is less than the last smallest distance (initially set as largest distance across the cube) and update if it is
                            if(distSqr <= minDistSqr_private){

                                NextIntersectFound = true;
                                minDistSqr_private = distSqr;
                                indClosest_private = j;

                            }

                        }
                    }

                    #pragma omp critical
                    {

                        if(minDistSqr_private <= minDistSqr){

                            minDistSqr = minDistSqr_private;
                            indClosest = indClosest_private;

                        }

                    }

                }

                // If found the next point, add the distance onto the total length of string, delete the intersection point being followed and adjust the next intersection point if needed

                if(NextIntersectFound){

                    if(splitLength){

                        if(xString[fInd] >= -0.25*nx*dx && xString[fInd] <= 0.25*nx*dx && yString[fInd] >= -0.25*ny*dy && yString[fInd] <= 0.25*ny*dy && zString[fInd] >= -0.25*nz*dz &&
                           zString[fInd] <= 0.25*nz*dz){ stringLength[0] += sqrt(minDistSqr); }
                        else{ stringLength[1] += sqrt(minDistSqr); }

                    } else{ stringLength[0] += sqrt(minDistSqr); }

                    xString.erase(xString.begin() + fInd);
                    yString.erase(yString.begin() + fInd);
                    zString.erase(zString.begin() + fInd);

                    if(indClosest<fInd){ fInd = indClosest; }
                    else{ fInd = indClosest - 1; } // if fInd is before indClosest, need to account for fInd being erased from the vector

                } else{

                    // Check if the reason there are no nearest neighbours is because a loop has been completed.

                    double distSqr = 0;

                    if(abs(xString[fInd] - beginCoords[0]) <= 0.5*nx*dx){ distSqr = pow(xString[fInd] - beginCoords[0],2); }
                    else if(xyBC == "periodic"){ distSqr = pow(nx*dx - abs(xString[fInd] - beginCoords[0]),2); }

                    if(abs(yString[fInd] - beginCoords[1]) <= 0.5*ny*dy){ distSqr += pow(yString[fInd] - beginCoords[1],2); }
                    else if(xyBC == "periodic"){ distSqr += pow(ny*dy - abs(yString[fInd] - beginCoords[1]),2); }

                    if(abs(zString[fInd] - beginCoords[2]) <= 0.5*nz*dz){ distSqr += pow(zString[fInd] - beginCoords[2],2); }
                    else if(zBC == "periodic"){ distSqr += pow(nz*dz - abs(zString[fInd] - beginCoords[2]),2); }

                    if(distSqr <= minDistSqr){ // If close enough to the starting point, count this length of string too

                        if(splitLength){

                            if(xString[fInd] >= -0.25*nx*dx && xString[fInd] <= 0.25*nx*dx && yString[fInd] >= -0.25*ny*dy && yString[fInd] <= 0.25*ny*dy && zString[fInd] >= -0.25*nz*dz &&
                               zString[fInd] <= 0.25*nz*dz){ stringLength[0] += sqrt(minDistSqr); }
                            else{ stringLength[1] += sqrt(minDistSqr); }

                        } else{ stringLength[0] += sqrt(minDistSqr); }

                    }

                    // Erase the point currently being following.

                    xString.erase(xString.begin() + fInd);
                    yString.erase(yString.begin() + fInd);
                    zString.erase(zString.begin() + fInd);

                    // Start the process again by resetting fInd to zero (which will now reference a new intersection point that hasn't been used yet). 

                    fInd = 0; // Reset fInd

                    // Store the coordinates of the new starting point
                    beginCoords[0] = xString[fInd];
                    beginCoords[1] = yString[fInd];
                    beginCoords[2] = zString[fInd];

                }

            }

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

        // If xyBC is absorbing, change loop to go over boundaries as this is how I am implementing the condition. Otherwise, don't change it from main loop.

        if(xyBC == "absorbing"){

            xEdge1 = 0;
            xEdge2 = nx;
            yEdge1 = 0;
            yEdge2 = ny;

        }

        #pragma omp parallel for default(none) shared(phi,theta,phitt,thetatt,tNow,tPast,xEdge1,xEdge2,yEdge1,yEdge2,zEdge1,zEdge2) private(j,k,comp)
        for(i=xEdge1;i<xEdge2;i++){
            for(j=yEdge1;j<yEdge2;j++){
                for(k=zEdge1;k<zEdge2;k++){

                    for(comp=0;comp<2;comp++){

                        phi(comp,tPast,i,j,k) = 2*phi(comp,tNow,i,j,k) - phi(comp,tPast,i,j,k) + dt*dt*phitt(comp,i,j,k);

                    }

                    for(comp=0;comp<3;comp++){

                        theta(comp,tPast,i,j,k) = 2*theta(comp,tNow,i,j,k) - theta(comp,tPast,i,j,k) + dt*dt*thetatt(comp,i,j,k);

                    }
                }
            }
        }

        if(xyBC == "periodic"){ deviationParameter = deviationParameter/(nx*ny); } // Include boundaries in average
        else{ deviationParameter = deviationParameter/((nx-2)*(ny-2)); }

        if(zBC == "periodic"){ deviationParameter = deviationParameter/nz; }
        else{ deviationParameter = deviationParameter/(nz-2); }


        valsPerLoop << energy << " " << deviationParameter;
        if(stringPos){ valsPerLoop << " " << pos; }
        if(stringDetect){ 

            valsPerLoop << " " << stringLength[0]; 
            if(splitLength){ valsPerLoop << " " << stringLength[1]; }

        }

        valsPerLoop << endl;


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

                finalField << phi(0,tPast,i,j,k) << " " << phi(1,tPast,i,j,k) << " " << theta(0,tPast,i,j,k) << " " << theta(1,tPast,i,j,k) << " " << theta(2,tPast,i,j,k) << endl;
                powerdensityOut << powerdensity(i,j,k) << endl;
                test1 << energydensity(tNow,i,j,k) << endl;

            }
        }
    }


    gettimeofday(&end,NULL);

    cout << "Time taken: " << end.tv_sec - start.tv_sec << "s" << endl;


	return 0;

}


