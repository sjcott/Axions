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
const int nt = 10001;
const double dx = 0.5;
const double dy = 0.5;
const double dz = 0.5;
const double dt = 0.05;

const double lambda = 1;
const double eta = 1;
const double g = 0.5;

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
	double phixx,phiyy,phizz,energy,deviationParameter,damp,phit[2],phiMagSqr,curx,cury,curz,Fxy_y,Fxz_z,Fyx_x,Fyz_z,Fzx_x,Fzy_y;
	int x0,y0,z0,i,j,k,TimeStep,gifStringPosFrame,tNow,tPast,counter,comp,imx,ipx,imy,ipy,imz,ipz,ipxmy,ipxmz,imxpy,ipymz,imxpz,imypz;
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

        if(TimeStep>counter and rank==0){

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

        	// Convert to global position in array to do modulo arithmetic. The second to last term gives ny*nz*floor((i+dataStart)/(ny*nz)). The last term converts back to the position in the local array 

        	imy = (i+dataStart-nz+ny*nz)%(ny*nz) + ( (i+dataStart)/(ny*nz) )*ny*nz - dataStart; 
        	ipy = (i+dataStart+nz)%(ny*nz) + ( (i+dataStart)/(ny*nz) )*ny*nz - dataStart;

        	imz = (i+dataStart-1+nz)%nz + ( (i+dataStart)/nz )*nz - dataStart;
        	ipz = (i+dataStart+1)%nz + ( (i+dataStart)/nz )*nz - dataStart;

        	// Additionally needed for wilson loop calculations. Avoid using x shifted points first as this makes the calculations more complicated and some of these points aren't in the correct positions

        	ipxmy = imy+ny*nz;
        	ipxmz = imz+ny*nz;
        	imxpy = ipy-ny*nz;
        	ipymz = (ipy+dataStart-1+nz)%nz + ( (ipy+dataStart)/nz )*nz - dataStart;
        	imxpz = ipz-ny*nz;
        	imypz = (imy+dataStart+1)%nz + ( (imy+dataStart)/nz )*nz - dataStart;

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


            // Calculate the derivatives of the field tensor (lattice version).

            Fxy_y = ( sin(theta[totSize*tNow+i] + theta[totSize*(nts+tNow)+ipx] - theta[totSize*tNow+ipy] - theta[totSize*(nts+tNow)+i]) - 
                      sin(theta[totSize*tNow+imy] + theta[totSize*(nts+tNow)+ipxmy] - theta[totSize*tNow+i] - theta[totSize*(nts+tNow)+imy]) )/(dy*dy);

            Fxz_z = ( sin(theta[totSize*tNow+i] + theta[totSize*(nts*2+tNow)+ipx] - theta[totSize*tNow+ipz] - theta[totSize*(nts*2+tNow)+i]) -
                      sin(theta[totSize*tNow+imz] + theta[totSize*(nts*2+tNow)+ipxmz] - theta[totSize*tNow+i] - theta[totSize*(nts*2+tNow)+imz]) )/(dz*dz);

            Fyx_x = ( sin(theta[totSize*(nts+tNow)+i] + theta[totSize*tNow+ipy] - theta[totSize*(nts+tNow)+ipx] - theta[totSize*tNow+i]) -
                      sin(theta[totSize*(nts+tNow)+imx] + theta[totSize*tNow+imxpy] - theta[totSize*(nts+tNow)+i] - theta[totSize*tNow+imx]) )/(dx*dx);

            Fyz_z = ( sin(theta[totSize*(nts+tNow)+i] + theta[totSize*(nts*2+tNow)+ipy] - theta[totSize*(nts+tNow)+ipz] - theta[totSize*(nts*2+tNow)+i]) -
                      sin(theta[totSize*(nts+tNow)+imz] + theta[totSize*(nts*2+tNow)+ipymz] - theta[totSize*(nts+tNow)+i] - theta[totSize*(nts*2+tNow)+imz]) )/(dz*dz);

            Fzx_x = ( sin(theta[totSize*(nts*2+tNow)+i] + theta[totSize*tNow+ipz] - theta[totSize*(nts*2+tNow)+ipx] - theta[totSize*tNow+i]) -
                      sin(theta[totSize*(nts*2+tNow)+imx] + theta[totSize*tNow+imxpz] - theta[totSize*(nts*2+tNow)+i] - theta[totSize*tNow+imx]) )/(dx*dx);

            Fzy_y = ( sin(theta[totSize*(nts*2+tNow)+i] + theta[totSize*(nts+tNow)+ipz] - theta[totSize*(nts*2+tNow)+ipy] - theta[totSize*(nts+tNow)+i]) -
                      sin(theta[totSize*(nts*2+tNow)+imy] + theta[totSize*(nts+tNow)+imypz] - theta[totSize*(nts*2+tNow)+i] - theta[totSize*(nts+tNow)+imy]) )/(dy*dy);

            thetatt[i-haloSize] = -2*g*g*curx - Fxy_y - Fxz_z - damp*( theta[totSize*tNow+i] - theta[totSize*tPast+i] )/dt;
            thetatt[coreSize+i-haloSize] = -2*g*g*cury - Fyx_x - Fyz_z - damp*( theta[totSize*(nts+tNow)+i] - theta[totSize*(nts+tPast)+i] )/dt;
            thetatt[coreSize*2+i-haloSize] = -2*g*g*curz - Fzx_x - Fzy_y - damp*( theta[totSize*(nts*2+tNow)+i] - theta[totSize*(nts*2+tPast)+i] )/dt;

        }

        // Update the core

        for(i=haloSize;i<coreSize+haloSize;i++){

        	for(comp=0;comp<2;comp++){

        		phi[totSize*(nts*comp+tPast)+i] = 2*phi[totSize*(nts*comp+tNow)+i] - phi[totSize*(nts*comp+tPast)+i] + dt*dt*phitt[coreSize*comp+i-haloSize];

        	}

        	for(comp=0;comp<3;comp++){

        		theta[totSize*(nts*comp+tPast)+i] = 2*theta[totSize*(nts*comp+tNow)+i] - theta[totSize*(nts*comp+tPast)+i] + dt*dt*thetatt[coreSize*comp+i-haloSize];

        	}

        }

        // Send sections of the core that are haloes for the other processes across to the relevant process. Then receive data for the halo of this process.

        for(comp=0;comp<2;comp++){ 

        	MPI_Send(&phi[totSize*(nts*comp+tPast)+haloSize],haloSize,MPI_DOUBLE,(rank-1+size)%size,0,MPI_COMM_WORLD);
        	MPI_Send(&phi[totSize*(nts*comp+tPast)+coreSize],haloSize,MPI_DOUBLE,(rank+1)%size,0,MPI_COMM_WORLD);

        	MPI_Recv(&phi[totSize*(nts*comp+tPast)],haloSize,MPI_DOUBLE,(rank-1+size)%size,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        	MPI_Recv(&phi[totSize*(nts*comp+tPast)+coreSize+haloSize],haloSize,MPI_DOUBLE,(rank+1)%size,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

        }

        for(comp=0;comp<3;comp++){ 

        	MPI_Send(&theta[totSize*(nts*comp+tPast)+haloSize],haloSize,MPI_DOUBLE,(rank-1+size)%size,0,MPI_COMM_WORLD);
        	MPI_Send(&theta[totSize*(nts*comp+tPast)+coreSize],haloSize,MPI_DOUBLE,(rank+1)%size,0,MPI_COMM_WORLD);

        	MPI_Recv(&theta[totSize*(nts*comp+tPast)],haloSize,MPI_DOUBLE,(rank-1+size)%size,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        	MPI_Recv(&theta[totSize*(nts*comp+tPast)+coreSize+haloSize],haloSize,MPI_DOUBLE,(rank+1)%size,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

        }

        // Code for testing the output

        if(TimeStep==nt-1){

	        if(rank==0){

	        	vector<double> phittOut(2*nPos,0.0), thetattOut(3*nPos,0.0);
	        	for(i=0;i<coreSize;i++){ phittOut[i] = phitt[i]; phittOut[nPos+i] = phitt[coreSize+i]; thetattOut[i] = thetatt[i]; thetattOut[nPos+i] = thetatt[coreSize+i]; thetattOut[2*nPos+i] = thetatt[2*coreSize+i]; }

	        	for(i=1;i<size;i++){

	        		int localCoreStart;
	        		int localCoreSize;
	        		if(i<chunkRem){ localCoreStart = i*(chunk+1); localCoreSize = chunk+1; }
	        		else{ localCoreStart = i*chunk + chunkRem; localCoreSize = chunk; }

	        		MPI_Recv(&phittOut[localCoreStart],localCoreSize,MPI_DOUBLE,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	        		MPI_Recv(&phittOut[nPos+localCoreStart],localCoreSize,MPI_DOUBLE,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

	        		MPI_Recv(&thetattOut[localCoreStart],localCoreSize,MPI_DOUBLE,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	        		MPI_Recv(&thetattOut[nPos+localCoreStart],localCoreSize,MPI_DOUBLE,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	        		MPI_Recv(&thetattOut[2*nPos+localCoreStart],localCoreSize,MPI_DOUBLE,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);  

				}

	        	for(i=0;i<nPos;i++){ testMerge << phittOut[i] << " " << phittOut[nPos+i] << " " << thetattOut[i] << " " << thetattOut[nPos+i] << " " << thetattOut[2*nPos+i] << endl; }

	        	// testArray[1] = 10;

	        	// cout << *foo << " " << *foo2 << endl;


	        } else{ MPI_Send(&phitt[0],coreSize,MPI_DOUBLE,0,0,MPI_COMM_WORLD); MPI_Send(&phitt[coreSize],coreSize,MPI_DOUBLE,0,0,MPI_COMM_WORLD); MPI_Send(&thetatt[0],coreSize,MPI_DOUBLE,0,0,MPI_COMM_WORLD); MPI_Send(&thetatt[coreSize],coreSize,MPI_DOUBLE,0,0,MPI_COMM_WORLD); MPI_Send(&thetatt[2*coreSize],coreSize,MPI_DOUBLE,0,0,MPI_COMM_WORLD); }

		}

    // Barrier before going to the next timestep. Not sure if strictly neccessary but I'm a paranoid man.

    MPI_Barrier(MPI_COMM_WORLD);

    }

    if(rank==0){

	    cout << "\rTimestep " << nt << " completed." << endl;

	    gettimeofday(&end,NULL);

	    cout << "Time taken: " << end.tv_sec - start.tv_sec << "s" << endl;

	}


    MPI_Finalize();


	return 0;

}