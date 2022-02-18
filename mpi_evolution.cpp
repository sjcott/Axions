#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <sys/time.h>
#include <complex>
#include <vector>
#include <mpi.h>

using namespace std;

// Never adjusted but useful to define for functions
const int nts = 2; // Number of time steps saved in data arrays

const int nx = 11;
const int ny = 11;
const int nz = 11;
const int nPos = nx*ny*nz;
const int nt = 1;
const double dx = 0.5;
const double dy = 0.5;
const double dz = 0.5;
const double dt = 0.05;

const double lambda = 1;
const double eta = 1;
const double g = 0;

const int damped_nt = 200; // Number of time steps for which damping is imposed. Useful for random initial conditions
const double dampFac = 0.5; // magnitude of damping term, unclear how strong to make this

// Below has not been implemented for gauge fields yet - so only works with global strings.
// Set alpha to zero to recover a non-expanding universe. Note that non-zero is not standard expansion but PRS algorithm.

const double alpha = 3; // Factor multiplying hubble damping term for use in PRS algorithm. alpha = #dims has been claimed to give similar dynamics without changing string width.
const double scaling = 1; // Power law scaling of the scale factor. Using conformal time so rad dom is gamma=1 while matter dom is gamma=2.

const bool makeGif = false; // Outputs data to make a gif of the isosurfaces
const bool makeStringPosGif = false; // Outputs data to make a gif of the calculated string position and curvature data
const int saveFreq = 5;
const int countRate = 10;

// Below has been removed and code now assumes periodic boundaries. Shouldn't be too tricky to add it back in if neccessary

//const string xyBC = "periodic"; // Allows for "neumann" (covariant derivatives set to zero), "absorbing", "periodic" or fixed (any other string will choose this option) boundary conditions.
//const string zBC = "periodic"; // Allows for "neumann" (covariant derivatives set to zero), "periodic" or fixed (any other string will choose this option) boundary conditions.

const bool stringPos = false;
const bool stationary_ic = true; // true if the initial conditions code doesn't also define the field values at second timestep.
const bool stringDetect = true; // true if you want the code to find which faces the strings pass through. May try to calculate string length later too.
const bool detectBuffer = true; // true if you want to don't want the code to find strings during damping process. Usually because random ic means you will find A LOT of strings --> memory issues.
const bool splitLength = false;
const bool finalOut = false;

int calcInd(int comp,int tStep,int pos){ return nPos*(nts*comp + tStep) + pos; }

int main(int argc, char ** argv){

	// Initialize MPI

    // Init MPI
    MPI_Init( &argc, &argv);

    // Get the rank and size
    int rank, size;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &size );

    if(size==1){cout << "Warning: Only one processor being used. This code is not designed for only one processor and may not work." << endl;}

    int chunk = nPos/size;
    int chunkRem = nPos - size*chunk;


    int coreSize;
	int haloSize = 2*nz*ny;
	if(rank>=chunkRem){ coreSize = chunk; }
	else{ coreSize = chunk+1; }
	int totSize = coreSize + 2*haloSize;

	cout << "Rank: " << rank << ", coreSize: " << coreSize << endl; 

	vector<double> phi(2*2*totSize, 0.0), theta(3*2*totSize, 0.0), phitt(2*coreSize, 0.0), thetatt(3*coreSize, 0.0), energydensity(coreSize, 0.0), gaussDeviation(coreSize, 0.0);
	double phixx,phiyy,phizz,energy,deviationParameter,damp,phit[2],phiMagSqr;
	int x0,y0,z0,i,j,k,TimeStep,gifStringPosFrame,tNow,tPast,counter,comp,imx,ipx,imy,ipy,imz,ipz;
	int c[2] = {1,-1}; // Useful definition to allow covariant deviative to be calculated when looping over components.

	struct timeval start, end;
	if(rank==0){ gettimeofday(&start, NULL); }

    string file_path = __FILE__;
    string dir_path = file_path.substr(0,file_path.find_last_of('/'));
    stringstream ss;

    MPI_Barrier(MPI_COMM_WORLD);

    string input;
    if(rank==0){

    	cout << "Enter a tag for output files: " << flush;
    	cin >> input;

    }

    MPI_Barrier(MPI_COMM_WORLD); // Allows all other processes to start once user input has been received.

    string icPath = dir_path + "/Data/ic.txt";
    string finalFieldPath = dir_path + "/Data/finalField.txt";
    string valsPerLoopPath = dir_path + "/Data/valsPerLoop_" + input + ".txt";
    string testPath = dir_path + "/Data/test_" + to_string(rank) + ".txt";
    string testMergePath = dir_path + "/Data/test.txt";

    ifstream ic (icPath.c_str());
    ofstream finalField (finalFieldPath.c_str());
    ofstream valsPerLoop (valsPerLoopPath.c_str());
    ofstream test (testPath.c_str());
    ofstream testMerge (testMergePath.c_str());

    x0 = int(0.5*(nx-1));
    y0 = int(0.5*(ny-1));
    z0 = int(0.5*(nz-1));

    int coreStart, coreEnd;
    if(rank < chunkRem){ coreStart = rank*(chunk+1); coreEnd = (rank+1)*(chunk+1); }
    else{ coreStart = rank*chunk+chunkRem; coreEnd = (rank+1)*chunk+chunkRem; }
    int dataStart = coreStart-haloSize;
    int dataEnd = coreEnd+haloSize;

    double wasteData[5];

    if(stationary_ic){

    	for(i=0;i<nPos;i++){

    		// Only assign it to the local array if this point belongs to the core or halo. Otherwise just waste it. Use modulus operator to deal with periodicity
    		// Uses index calculation totSize*(nTimeSteps*component + TimeStep) + pos

    		if(dataStart<0){ // Need to deal with periodicity at the start of the data array

    			if(i>=(dataStart+nPos)%nPos){ // Backward halo

    				int arrayPos = i-(dataStart+nPos)%nPos;

    				ic >> phi[arrayPos] >> phi[totSize*nts+arrayPos] >> theta[arrayPos] >> theta[totSize*nts+arrayPos] >> theta[totSize*nts*2+arrayPos];

    				// Second time step is equal to the first

    				phi[totSize+arrayPos] = phi[arrayPos];
    				phi[totSize*(nts+1)+arrayPos] = phi[totSize*nts+arrayPos];
    				theta[totSize+arrayPos] = theta[arrayPos];
    				theta[totSize*(nts+1)+arrayPos] = theta[totSize*nts+arrayPos];
    				theta[totSize*(nts*2+1)+arrayPos] = theta[totSize*nts*2+arrayPos];

    			} else if(i<dataEnd){ // The rest of the data

    				int arrayPos = i+haloSize; // Shift across to account for the backward halo at the start

    				ic >> phi[arrayPos] >> phi[totSize*nts+arrayPos] >> theta[arrayPos] >> theta[totSize*nts+arrayPos] >> theta[totSize*nts*2+arrayPos];

    				phi[totSize+arrayPos] = phi[arrayPos];
    				phi[totSize*(nts+1)+arrayPos] = phi[totSize*nts+arrayPos];
    				theta[totSize+arrayPos] = theta[arrayPos];
    				theta[totSize*(nts+1)+arrayPos] = theta[totSize*nts+arrayPos];
    				theta[totSize*(nts*2+1)+arrayPos] = theta[totSize*nts*2+arrayPos];

    			} else{ ic >> wasteData[0] >> wasteData[1] >> wasteData[2] >> wasteData[3] >> wasteData[4]; } // Don't need these so waste them into an unused variable 


    		} else if(dataEnd>nPos){ // Need to deal with periodicity at the end of the data array

    			if(i>=dataStart){ // All of the array except for the forward halo

    				int arrayPos = i-dataStart;

					ic >> phi[arrayPos] >> phi[totSize*nts+arrayPos] >> theta[arrayPos] >> theta[totSize*nts+arrayPos] >> theta[totSize*nts*2+arrayPos];

    				phi[totSize+arrayPos] = phi[arrayPos];
    				phi[totSize*(nts+1)+arrayPos] = phi[totSize*nts+arrayPos];
    				theta[totSize+arrayPos] = theta[arrayPos];
    				theta[totSize*(nts+1)+arrayPos] = theta[totSize*nts+arrayPos];
    				theta[totSize*(nts*2+1)+arrayPos] = theta[totSize*nts*2+arrayPos]; 				

    			} else if(i<dataEnd%nPos){ // The forward halo

    				int arrayPos = i+coreSize+haloSize;

    				ic >> phi[arrayPos] >> phi[totSize*nts+arrayPos] >> theta[arrayPos] >> theta[totSize*nts+arrayPos] >> theta[totSize*nts*2+arrayPos];

    				phi[totSize+arrayPos] = phi[arrayPos];
    				phi[totSize*(nts+1)+arrayPos] = phi[totSize*nts+arrayPos];
    				theta[totSize+arrayPos] = theta[arrayPos];
    				theta[totSize*(nts+1)+arrayPos] = theta[totSize*nts+arrayPos];
    				theta[totSize*(nts*2+1)+arrayPos] = theta[totSize*nts*2+arrayPos]; 

    			} else{ ic >> wasteData[0] >> wasteData[1] >> wasteData[2] >> wasteData[3] >> wasteData[4]; }

    		} else{ // In the middle of the array so don't need to deal with periodicity

    			if(i>=dataStart and i<dataEnd){

    				int arrayPos = i-dataStart;

    				ic >> phi[arrayPos] >> phi[totSize*nts+arrayPos] >> theta[arrayPos] >> theta[totSize*nts+arrayPos] >> theta[totSize*nts*2+arrayPos];

    				phi[totSize+arrayPos] = phi[arrayPos];
    				phi[totSize*(nts+1)+arrayPos] = phi[totSize*nts+arrayPos];
    				theta[totSize+arrayPos] = theta[arrayPos];
    				theta[totSize*(nts+1)+arrayPos] = theta[totSize*nts+arrayPos];
    				theta[totSize*(nts*2+1)+arrayPos] = theta[totSize*nts*2+arrayPos]; 

    			} else{ ic >> wasteData[0] >> wasteData[1] >> wasteData[2] >> wasteData[3] >> wasteData[4]; }

    		}

    	}

    } else{

    	for(TimeStep=0;TimeStep<2;TimeStep++){
    		for(i=0;i<nPos;i++){

    			if(dataStart<0){ // Need to deal with periodicity at the start of the data array

	    			if(i>=(dataStart+nPos)%nPos){ // Backward halo

	    				int arrayPos = i-(dataStart+nPos)%nPos;

	    				ic >> phi[totSize*TimeStep+arrayPos] >> phi[totSize*(nts+TimeStep)+arrayPos] >> theta[totSize*TimeStep+arrayPos] >> theta[totSize*(nts+TimeStep)+arrayPos] 
	    				   >> theta[totSize*(nts*2+TimeStep)+arrayPos];

	    			} else if(i<dataEnd){ // The rest of the data

	    				int arrayPos = i+haloSize; // Shift across to account for the backward halo at the start

	    				ic >> phi[totSize*TimeStep+arrayPos] >> phi[totSize*(nts+TimeStep)+arrayPos] >> theta[totSize*TimeStep+arrayPos] >> theta[totSize*(nts+TimeStep)+arrayPos] 
	    				   >> theta[totSize*(nts*2+TimeStep)+arrayPos];

	    			} else{ ic >> wasteData[0] >> wasteData[1] >> wasteData[2] >> wasteData[3] >> wasteData[4]; } // Don't need these so waste them into an unused variable 


	    		} else if(dataEnd>nPos){ // Need to deal with periodicity at the end of the data array

	    			if(i>=dataStart){ // All of the array except for the forward halo

	    				int arrayPos = i-dataStart;

						ic >> phi[totSize*TimeStep+arrayPos] >> phi[totSize*(nts+TimeStep)+arrayPos] >> theta[totSize*TimeStep+arrayPos] >> theta[totSize*(nts+TimeStep)+arrayPos] 
	    				   >> theta[totSize*(nts*2+TimeStep)+arrayPos];				

	    			} else if(i<dataEnd%nPos){ // The forward halo

	    				int arrayPos = i+coreSize+haloSize;

	    				ic >> phi[totSize*TimeStep+arrayPos] >> phi[totSize*(nts+TimeStep)+arrayPos] >> theta[totSize*TimeStep+arrayPos] >> theta[totSize*(nts+TimeStep)+arrayPos] 
	    				   >> theta[totSize*(nts*2+TimeStep)+arrayPos];

	    			} else{ ic >> wasteData[0] >> wasteData[1] >> wasteData[2] >> wasteData[3] >> wasteData[4]; }

	    		} else{ // In the middle of the array so don't need to deal with periodicity

	    			if(i>=dataStart and i<dataEnd){

	    				int arrayPos = i-dataStart;

	    				ic >> phi[totSize*TimeStep+arrayPos] >> phi[totSize*(nts+TimeStep)+arrayPos] >> theta[totSize*TimeStep+arrayPos] >> theta[totSize*(nts+TimeStep)+arrayPos] 
	    				   >> theta[totSize*(nts*2+TimeStep)+arrayPos];

	    			} else{ ic >> wasteData[0] >> wasteData[1] >> wasteData[2] >> wasteData[3] >> wasteData[4]; }

	    		}

    		}
    	}

    }

    // All relevant data loaded.

    gifStringPosFrame = 0;
    counter = 0;

    for(TimeStep=0;TimeStep<nt;TimeStep++){

        double time = 1 + TimeStep*dt; // Conformal time, starting at eta = 1.

        if(TimeStep>counter){

            cout << "\rTimestep " << TimeStep-1 << " completed." << flush;

            counter += countRate;

        }

        // Is damping switched on or not?
        if(TimeStep<damped_nt){ damp = dampFac; }
        else{ damp = 0; }

        tNow = (TimeStep+1)%2;
        tPast = TimeStep%2;


        // Calculate time derivatives using EoMs

        energy = 0;
        deviationParameter = 0;

        for(i=haloSize;i<coreSize+haloSize;i++){ // Now evolve the core data

        	// No need to worry about periodicity with the x neighbours because halo is designed to contain them

        	imx = i-ny*nz;
        	ipx = i+ny*nz;

        	// Need to account for the periodicity of the space for the other two directions

        	imy = (i+dataStart-nz+ny*nz)%(ny*nz) + ( (i+dataStart)/(ny*nz) )*ny*nz - dataStart; // Last term gives ny*nz*floor(i/(ny*nz))
        	ipy = (i+dataStart+nz)%(ny*nz) + ( (i+dataStart)/(ny*nz) )*ny*nz - dataStart;

        	imz = (i+dataStart-1+nz)%nz + ( (i+dataStart)/nz )*nz - dataStart;
        	ipz = (i+dataStart+1)%nz + ( (i+dataStart)/nz )*nz - dataStart;

	        phiMagSqr = pow(phi[totSize*tNow+i],2) + pow(phi[totSize*(nts+tNow)+i],2);

            // Loop over phi components

            for(comp=0;comp<2;comp++){

                 // 2nd order spatial derivatives calculated with 2nd order finite difference

                // c is 1 when comp = 0 and -1 when comp = 1

                phixx = ( cos(theta[totSize*tNow+i])*phi[totSize*(nts*comp+tNow)+ipx] + c[comp]*sin(theta[totSize*tNow+i])*phi[totSize*(nts*(comp+c[comp])+tNow)+ipx] - 2*phi[totSize*(nts*comp+tNow)+i]
                		+ cos(theta[totSize*tNow+imx])*phi[totSize*(nts*comp+tNow)+imx] - c[comp]*sin(theta[totSize*tNow+imx])*phi[totSize*(nts*(comp+c[comp])+tNow)+imx] )/(dx*dx);

                phiyy = ( cos(theta[totSize*(nts+tNow)+i])*phi[totSize*(nts*comp+tNow)+ipy] + c[comp]*sin(theta[totSize*(nts+tNow)+i])*phi[totSize*(nts*(comp+c[comp])+tNow)+ipy] - 2*phi[totSize*(nts*comp+tNow)+i]
                		+ cos(theta[totSize*(nts+tNow)+imy])*phi[totSize*(nts*comp+tNow)+imy] - c[comp]*sin(theta[totSize*(nts+tNow)+imy])*phi[totSize*(nts*(comp+c[comp])+tNow)+imy] )/(dy*dy);

                phizz = ( cos(theta[totSize*(nts*2+tNow)+i])*phi[totSize*(nts*comp+tNow)+ipz] + c[comp]*sin(theta[totSize*(nts*2+tNow)+i])*phi[totSize*(nts*(comp+c[comp])+tNow)+ipz] - 2*phi[totSize*(nts*comp+tNow)+i]
                		+ cos(theta[totSize*(nts*2+tNow)+imz])*phi[totSize*(nts*comp+tNow)+imz] - c[comp]*sin(theta[totSize*(nts*2+tNow)+imz])*phi[totSize*(nts*(comp+c[comp])+tNow)+imz] )/(dz*dz);


                phit[comp] = ( phi[totSize*(nts*comp+tNow)+i] - phi[totSize*(nts*comp+tPast)+i] )/dt;

                // Calculate second order time derivatives

                phitt[coreSize*comp+i-haloSize] = phixx + phiyy + phizz - 0.5*lambda*(phiMagSqr - pow(eta,2))*phi[totSize*(nts*comp+tNow)+i] - (damp + alpha*scaling/time)*phit[comp];


                // Calculate first order derivatives for energy

                // phix[comp] = ( cos(theta(0,tNow,i,j,k))*phi(comp,tNow,i+1,j,k) + c[comp]*sin(theta(0,tNow,i,j,k))*phi(comp+c[comp],tNow,i+1,j,k) 
                //              - cos(theta(0,tNow,i-1,j,k))*phi(comp,tNow,i-1,j,k) + c[comp]*sin(theta(0,tNow,i-1,j,k))*phi(comp+c[comp],tNow,i-1,j,k) )/(2*dx);
                // phiy[comp] = ( cos(theta(1,tNow,i,j,k))*phi(comp,tNow,i,j+1,k) + c[comp]*sin(theta(1,tNow,i,j,k))*phi(comp+c[comp],tNow,i,j+1,k) 
                //              - cos(theta(1,tNow,i,j-1,k))*phi(comp,tNow,i,j-1,k) + c[comp]*sin(theta(1,tNow,i,j-1,k))*phi(comp+c[comp],tNow,i,j-1,k) )/(2*dy);

                //Try without central differencing

                // phix[comp] = ( phi(comp,tNow,i,j,k) - cos(theta(0,tNow,im,j,k))*phi(comp,tNow,im,j,k) + c[comp]*sin(theta(0,tNow,im,j,k))*phi(comp+c[comp],tNow,im,j,k) )/dx;
                // phiy[comp] = ( phi(comp,tNow,i,j,k) - cos(theta(1,tNow,i,jm,k))*phi(comp,tNow,i,jm,k) + c[comp]*sin(theta(1,tNow,i,jm,k))*phi(comp+c[comp],tNow,i,jm,k) )/dy;
                // phiz[comp] = ( phi(comp,tNow,i,j,k) - cos(theta(2,tNow,i,j,km))*phi(comp,tNow,i,j,km) + c[comp]*sin(theta(2,tNow,i,j,km))*phi(comp+c[comp],tNow,i,j,km) )/dz;


            }


            // Calculate the currents

            curx = cos(theta[totSize*tNow+i])*(phi[totSize*(nts+tNow)+i]*phi[totSize*tNow+ipx] - phi[totSize*tNow+i]*phi[totSize*(nts+tNow)+ipx]) +
                   sin(theta[totSize*tNow+i])*(phi[totSize*tNow+i]*phi[totSize*tNow+ipx] + phi[totSize*(nts+tNow)+i]*phi[totSize*(nts+tNow)+ipx]);

            cury = cos(theta[totSize*(nts+tNow)+i])*(phi[totSize*(nts+tNow)+i]*phi[totSize*tNow+ipy] - phi[totSize*tNow+i]*phi[totSize*(nts+tNow)+ipy]) +
                   sin(theta[totSize*(nts+tNow)+i])*(phi[totSize*tNow+i]*phi[totSize*tNow+ipy] + phi[totSize*(nts+tNow)+i]*phi[totSize*(nts+tNow)+ipy]);

            curz = cos(theta[totSize*(nts*2+tNow)+i])*(phi[totSize*(nts+tNow)+i]*phi[totSize*tNow+ipz] - phi[totSize*tNow+i]*phi[totSize*(nts+tNow)+ipz]) +
                   sin(theta[totSize*(nts*2+tNow)+i])*(phi[totSize*tNow+i]*phi[totSize*tNow+ipz] + phi[totSize*(nts+tNow)+i]*phi[totSize*(nts+tNow)+ipz]);


            // Calculate the derivatives of the field tensor (lattice version). Additionally need ipxmy,ipxmz,imxpy

            Fxy_y = ( sin(theta[totSize*tNow+i] + theta[totSize*(nts+tNow)+ipx] - theta[totSize*tNow+ipy] - theta[totSize*(nts+tNow)+i]) - 
                      sin(theta[totSize*tNow+imy] + theta[totSize*(nts+tNow)+ipxmy] - theta[totSize*tNow+i] - theta[totSize*(nts+tNow)+imy]) )/(dy*dy);

            Fxz_z = ( sin(theta[totSize*tNow+i] + theta[totSize*(nts*2+tNow)+ipx] - theta[totSize*tNow+ipz] - theta[totSize*(nts*2+tNow)+i]) -
                      sin(theta[totSize*tNow+imz] + theta[totSize*(nts*2+tNow)+ipxmz] - theta[totSize*tNow+i] - theta[totSize*(nts*2+tNow)+imz]) )/(dz*dz);

            Fyx_x = ( sin(theta[totSize*(nts+tNow)+i] + theta[totSize*tNow+ipy] - theta[totSize*(nts+tNow)+ipx] - theta[totSize*tNow+i]) -
                      sin(theta[totSize*(nts+tNow)+imx] + theta[totSize*tNow+imxpy] - theta[totSize*(nts+tNow)+i] - theta[totSize*tNow+imx]) )/(dx*dx);

            Fyz_z = ( sin(theta[totSize*(nts+tNow)+i] + theta[totSize*(nts*2+tNow)+ipy] - theta[totSize*(nts+tNow)+ipz] - theta(2,tNow,i,j,k)) -
                      sin(theta(1,tNow,i,j,km) + theta(2,tNow,i,jp,km) - theta(1,tNow,i,j,k) - theta(2,tNow,i,j,km)) )/(dz*dz);

            Fzx_x = ( sin(theta(2,tNow,i,j,k) + theta(0,tNow,i,j,kp) - theta(2,tNow,ip,j,k) - theta(0,tNow,i,j,k)) -
                      sin(theta(2,tNow,im,j,k) + theta(0,tNow,im,j,kp) - theta(2,tNow,i,j,k) - theta(0,tNow,im,j,k)) )/(dx*dx);

            Fzy_y = ( sin(theta(2,tNow,i,j,k) + theta(1,tNow,i,j,kp) - theta(2,tNow,i,jp,k) - theta(1,tNow,i,j,k)) -
                      sin(theta(2,tNow,i,jm,k) + theta(1,tNow,i,jm,kp) - theta(2,tNow,i,j,k) - theta(1,tNow,i,jm,k)) )/(dy*dy);

            thetatt(0,i,j,k) = -2*g*g*curx - Fxy_y - Fxz_z - damp*( theta(0,tNow,i,j,k) - theta(0,tPast,i,j,k) )/dt;
            thetatt(1,i,j,k) = -2*g*g*cury - Fyx_x - Fyz_z - damp*( theta(1,tNow,i,j,k) - theta(1,tPast,i,j,k) )/dt;
            thetatt(2,i,j,k) = -2*g*g*curz - Fzx_x - Fzy_y - damp*( theta(2,tNow,i,j,k) - theta(2,tPast,i,j,k) )/dt;

        }

        if(rank==0){

        	vector<double> phixxOut(2*nPos,0.0), phiyyOut(2*nPos,0.0), phizzOut(2*nPos,0.0), phittOut(2*nPos,0.0);
        	for(i=0;i<coreSize;i++){ phixxOut[i] = phixxVec[i]; phixxOut[nPos+i] = phixxVec[coreSize+i]; phiyyOut[i] = phiyyVec[i]; phiyyOut[nPos+i] = phiyyVec[coreSize+i]; phizzOut[i] = phizzVec[i]; phizzOut[nPos+i] = phizzVec[coreSize+i]; phittOut[i] = phitt[i]; phittOut[nPos+i] = phitt[coreSize+i]; }

        	for(i=1;i<size;i++){

        		int localCoreStart;
        		int localCoreSize;
        		if(i<chunkRem){ localCoreStart = i*(chunk+1); localCoreSize = chunk+1; }
        		else{ localCoreStart = i*chunk + chunkRem; localCoreSize = chunk; }

				MPI_Recv(&phixxOut[localCoreStart],localCoreSize,MPI_DOUBLE,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				MPI_Recv(&phixxOut[nPos+localCoreStart],localCoreSize,MPI_DOUBLE,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE); 

				MPI_Recv(&phiyyOut[localCoreStart],localCoreSize,MPI_DOUBLE,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				MPI_Recv(&phiyyOut[nPos+localCoreStart],localCoreSize,MPI_DOUBLE,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE); 

				MPI_Recv(&phizzOut[localCoreStart],localCoreSize,MPI_DOUBLE,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				MPI_Recv(&phizzOut[nPos+localCoreStart],localCoreSize,MPI_DOUBLE,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

				MPI_Recv(&phittOut[localCoreStart],localCoreSize,MPI_DOUBLE,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				MPI_Recv(&phittOut[nPos+localCoreStart],localCoreSize,MPI_DOUBLE,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);  

			}

        	for(i=0;i<nPos;i++){ testMerge << phixxOut[i] << " " << phixxOut[nPos+i] << " " << phiyyOut[i] << " " << phiyyOut[nPos+i] << " " << phizzOut[i] << " " << phizzOut[nPos+i] << " " << phittOut[i] << " " << phittOut[nPos+i] << endl; }

        	// testArray[1] = 10;

        	// cout << *foo << " " << *foo2 << endl;


        } else{ MPI_Send(&phixxVec[0],coreSize,MPI_DOUBLE,0,0,MPI_COMM_WORLD); MPI_Send(&phixxVec[coreSize],coreSize,MPI_DOUBLE,0,0,MPI_COMM_WORLD); MPI_Send(&phiyyVec[0],coreSize,MPI_DOUBLE,0,0,MPI_COMM_WORLD); MPI_Send(&phiyyVec[coreSize],coreSize,MPI_DOUBLE,0,0,MPI_COMM_WORLD); MPI_Send(&phizzVec[0],coreSize,MPI_DOUBLE,0,0,MPI_COMM_WORLD); MPI_Send(&phizzVec[coreSize],coreSize,MPI_DOUBLE,0,0,MPI_COMM_WORLD); MPI_Send(&phitt[0],coreSize,MPI_DOUBLE,0,0,MPI_COMM_WORLD); MPI_Send(&phitt[coreSize],coreSize,MPI_DOUBLE,0,0,MPI_COMM_WORLD); }


        //MPI_Send( message_r1, 13, MPI_CHAR, 1, 0, MPI_COMM_WORLD );


    }


 //    if(rank==0){

 //    	int coreSize;
 //    	int haloSize = 2*nz*ny;
 //    	if(rank>=chunkRem){ coreSize = chunk; }
 //    	else{ coreSize = chunk+1; }
 //    	int totSize = coreSize + haloSize;

 //    	cout << "Rank: " << rank << ", coreSize: " << coreSize << endl; 

	// 	vector<double> phi(2*2*totSize, 0.0), theta(3*2*totSize, 0.0), phitt(2*totSize, 0.0), thetatt(3*totSize, 0.0), energydensity(totSize, 0.0), gaussDeviation(totSize, 0.0);
	// 	double energy,deviationParameter,damp;
	// 	int x0,y0,z0,i,j,k,TimeStep,gifStringPosFrame,tNow,tPast,counter;

	// 	struct timeval start, end;
	//     gettimeofday(&start, NULL);

	//     string file_path = __FILE__;
	//     string dir_path = file_path.substr(0,file_path.find_last_of('/'));
	//     stringstream ss;

	//     string input;
	//     cout << "Enter a tag for output files: " << flush;
	//     cin >> input;

	//     MPI_Barrier(MPI_COMM_WORLD); // Allows all other processes to start once user input has been received.

	//     string icPath = dir_path + "/Data/ic.txt";
	//     string finalFieldPath = dir_path + "/Data/finalField.txt";
	//     string valsPerLoopPath = dir_path + "/Data/valsPerLoop_" + input + ".txt";
	//     string testPath = dir_path + "/Data/test.txt";

	//     ifstream ic (icPath.c_str());
	//     ofstream finalField (finalFieldPath.c_str());
	//     ofstream valsPerLoop (valsPerLoopPath.c_str());
	//     ofstream test (testPath.c_str());

	//     x0 = int(0.5*(nx-1));
	//     y0 = int(0.5*(ny-1));
	//     z0 = int(0.5*(nz-1));

	//     double remCont1, remCont2;
	//     if(rank < chunkRem){ remCont1 = rank; remCont2 = rank+1; }
	//     else{ remCont1 = chunkRem; remCont2 = chunkRem; }
	//     double wasteData[5];

	//     if(stationary_ic){

	//     	for(i=0;i<nx*ny*nz;i++){

	//     		// Only assign it to the local array if this point belongs to the core or halo. Otherwise just waste it.

	//     		// Master process always starts at 0 and has to deal with the halo being wrapped to the other side of the array

	//     		if(i<(rank+1)*chunk+remCont2 + ny*nz){

	//     			// Calculating index with totSize*(nTimeSteps*comp + TimeStep) + pos and shifting across by the size of the first half of the halo

	//     			ic >> phi[i+ny*nz] >> phi[totSize*2*1+i+ny*nz] >> theta[i+ny*nz] >> theta[totSize*2*1+i+ny*nz] >> theta[totSize*2*2+i+ny*nz];

	//     			// Second time step is equal to the first

	//     			phi[totSize+i+ny*nz] = phi[i+ny*nz];
	//     			phi[totSize*(2*1+1)+i+ny*nz] = phi[totSize*2*1+i+ny*nz];
	//     			theta[totSize+i+ny*nz] = theta[i+ny*nz];
	//     			theta[totSize*(2*1+1)+i+ny*nz] = theta[totSize*2*1+i+ny*nz];
	//     			theta[totSize*(2*2+1)+i+ny*nz] = theta[totSize*2*2+i+ny*nz];

	//     			//test << i+ny*nz << endl << totSize*2*1+i+ny*nz << endl << totSize*2*2+i+ny*nz << endl << totSize+i+ny*nz << endl << totSize*(2*1+1)+i+ny*nz << endl << totSize*(2*2+1)+i+ny*nz << endl;

	//     		} else if(i>=nx*ny*nz - ny*nz){

	//     			// Same index calculation as above but want to now allocate the first half of the halo to the beginning of the array

	//     			ic >> phi[i-nx*ny*nz+ny*nz] >> phi[totSize*2*1+i-nx*ny*nz+ny*nz] >> theta[i-nx*ny*nz+ny*nz] >> theta[totSize*2*1+i-nx*ny*nz+ny*nz] >> theta[totSize*2*2+i-nx*ny*nz+ny*nz];

	//     			phi[totSize+i-nx*ny*nz+ny*nz] = phi[i-nx*ny*nz+ny*nz];
	//     			phi[totSize*(2*1+1)+i-nx*ny*nz+ny*nz] = phi[totSize*2*1+i-nx*ny*nz+ny*nz];
	//     			theta[totSize+i+ny*nz] = theta[i+ny*nz];
	//     			theta[totSize*(2*1+1)+i-nx*ny*nz+ny*nz] = theta[totSize*2*1+i-nx*ny*nz+ny*nz];
	//     			theta[totSize*(2*2+1)+i-nx*ny*nz+ny*nz] = theta[totSize*2*2+i-nx*ny*nz+ny*nz];

	//     			//test << i-nx*ny*nz+ny*nz << endl << totSize*2*1+i-nx*ny*nz+ny*nz << endl << totSize*2*2+i-nx*ny*nz+ny*nz << endl << totSize+i-nx*ny*nz+ny*nz << endl << totSize*(2*1+1)+i-nx*ny*nz+ny*nz << endl << totSize*(2*2+1)+i-nx*ny*nz+ny*nz << endl;


	//     		} else{

	//     			// Don't need these points so just throw them away into an unused variable

	//     			ic >> wasteData[0] >> wasteData[1] >> wasteData[2] >> wasteData[3] >> wasteData[4];

	//     		}

	//     	}

	//     } else{

	//     	for(TimeStep=0;TimeStep<2;TimeStep++){
	//     		for(i=0;i<nx*ny*nz;i++){

	//     			if(i<(rank+1)*chunk+remCont2 + ny*nz){

	//     				ic >> phi[totSize*TimeStep+i+ny*nz] >> phi[totSize*(2*1+TimeStep)+i+ny*nz] >> theta[totSize*TimeStep+i+ny*nz] >> theta[totSize*(2*1+TimeStep)+i+ny*nz] >> theta[totSize*(2*2+TimeStep)+i+ny*nz];

	//     			} else if(i>=nx*ny*nz - ny*nz){

	//     				ic >> phi[totSize*TimeStep+i-nx*ny*nz+ny*nz] >> phi[totSize*(2*1+TimeStep)+i-nx*ny*nz+ny*nz] >> theta[totSize*TimeStep+i-nx*ny*nz+ny*nz]
	//     				   >> theta[totSize*(2*1+TimeStep)+i-nx*ny*nz+ny*nz] >> theta[totSize*(2*2+TimeStep)+i-nx*ny*nz+ny*nz];

	//     			}

	//     		}
	//     	}

	//     }

	//     // All relevant data loaded. Now evolve these points.

	//     gifStringPosFrame = 0;
	//     counter = 0;

	//     for(TimeStep=0;TimeStep<1;TimeStep++){

	//         double time = 1 + TimeStep*dt; // Conformal time, starting at eta = 1.

	//         if(TimeStep>counter){

	//             cout << "\rTimestep " << TimeStep-1 << " completed." << flush;

	//             counter += countRate;

	//         }

	//         // Is damping switched on or not?
	//         if(TimeStep<damped_nt){ damp = dampFac; }
	//         else{ damp = 0; }

	//         tNow = (TimeStep+1)%2;
	//         tPast = TimeStep%2;


	//         // Calculate time derivatives using EoMs

	//         energy = 0;
	//         deviationParameter = 0;


	//         //MPI_Send( message_r1, 13, MPI_CHAR, 1, 0, MPI_COMM_WORLD );


	//     }

	// } else{

	// 	int coreSize;
 //    	int haloSize = 2*nz*ny;
 //    	if(rank>=chunkRem){ coreSize = chunk; }
 //    	else{ coreSize = chunk+1; }
 //    	int totSize = coreSize + haloSize;

 //    	cout << "Rank: " << rank << ", coreSize: " << coreSize << endl; 

	// 	vector<double> phi(2*2*totSize, 0.0), theta(3*2*totSize, 0.0), phitt(2*totSize, 0.0), thetatt(3*totSize, 0.0), energydensity(totSize, 0.0), gaussDeviation(totSize, 0.0);

	// 	MPI_Barrier(MPI_COMM_WORLD); // Wait for user input from master process.

	// 	string file_path = __FILE__;
	// 	string dir_path = file_path.substr(0,file_path.find_last_of('/'));
	// 	string icPath = dir_path + "/Data/ic.txt";
	// 	ifstream ic (icPath.c_str());

	// 	int i,j,k,TimeStep;

	//     if(stationary_ic){

	//         for(i=0;i<nx;i++){
	//             for(j=0;j<ny;j++){
	//                 for(k=0;k<nz;k++){

	//                 	double inData[5];

	//                     ic >> inData[0] >> inData[1] >> inData[2] >> inData[3] >> inData[4];

	//                 }
	//             }
	//         }

	//     } else{

	//         for(TimeStep=0;TimeStep<2;TimeStep++){
	//             for(i=0;i<nx;i++){
	//                 for(j=0;j<ny;j++){
	//                     for(k=0;k<nz;k++){

	//                         //ic >> phi[calcInd(0,TimeStep,i,j,k)] >> phi[calcInd(1,TimeStep,i,j,k)] >> theta[calcInd(0,TimeStep,i,j,k)] >> theta[calcInd(1,TimeStep,i,j,k)] >> theta[calcInd(2,TimeStep,i,j,k)];

	//                     }
	//                 }
	//             }
	//         }

	//     }


	// }

    MPI_Finalize();


	return 0;

}