#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <sys/time.h>
#include <vector>
#include <mpi.h>
#include <random>

using namespace std;

// Never adjusted but useful to define for functions
const int nts = 2; // Number of time steps saved in data arrays

const int nx = 401;
const int ny = 401;
const int nz = 401;
const int nPos = nx*ny*nz;
const int nt = 8001;
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

const bool makeGif = false; // Outputs data to make a gif of the isosurfaces. Not implemented in this version yet.
const bool makeStringPosGif = true; // Outputs data to make a gif of the calculated string position and curvature data
const int saveFreq = 5;
const int countRate = 10;
const string gifTag = "global_PRS_dx0p5_";

// How are the initial conditions generated? Below are all parameters used for generating (or loading) the initial conditions
const string ic_type = "random"; // Current options are data, stationary data and random
const double seed = 42;
const double mean = 0; // Probably always want this to be zero
const double stdev = 0.5;

// Below has been removed and code now assumes periodic boundaries. Shouldn't be too tricky to add it back in if neccessary

//const string xyBC = "periodic"; // Allows for "neumann" (covariant derivatives set to zero), "absorbing", "periodic" or fixed (any other string will choose this option) boundary conditions.
//const string zBC = "periodic"; // Allows for "neumann" (covariant derivatives set to zero), "periodic" or fixed (any other string will choose this option) boundary conditions.

const bool stringPos = false;
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

    int chunk = nPos/size;
    int chunkRem = nPos - size*chunk;

    int coreSize;
	int haloSize = 2*nz*ny;
	if(rank>=chunkRem){ coreSize = chunk; }
	else{ coreSize = chunk+1; }
	int totSize = coreSize + 2*haloSize;

    // Warnings

    if(rank==0){

    	if(size==1){ cout << "Warning: Only one processor being used. This code is not designed for only one processor and may not work." << endl; }
    	if(chunk<haloSize){ cout << "Warning: Chunk size is less than the halo size (i.e chunk neighbour data). Code currently assumes this is not the case so it probably won't work." << endl; }

    }

	vector<double> phi(2*2*totSize, 0.0), theta(3*2*totSize, 0.0), phitt(2*coreSize, 0.0), thetatt(3*coreSize, 0.0), energydensity(coreSize, 0.0), gaussDeviation(coreSize, 0.0);
	double phixx,phiyy,phizz,energy,deviationParameter,damp,phit[2],phiMagSqr,curx,cury,curz,Fxy_y,Fxz_z,Fyx_x,Fyz_z,Fzx_x,Fzy_y,localSeed;
	int x0,y0,z0,i,j,k,TimeStep,gifStringPosFrame,tNow,tPast,counter,comp,imx,ipx,imy,ipy,imz,ipz,ipxmy,ipxmz,imxpy,ipymz,imxpz,imypz;
	int c[2] = {1,-1}; // Useful definition to allow covariant deviative to be calculated when looping over components.

	struct timeval start, end;
	if(rank==0){ gettimeofday(&start, NULL); }

    string file_path = __FILE__;
    string dir_path = file_path.substr(0,file_path.find_last_of('/'));
    stringstream ss;

    MPI_Barrier(MPI_COMM_WORLD);

    // string input;
    // if(rank==0){

    // 	cout << "Enter a tag for output files: " << flush;
    // 	cin >> input;

    // }

    MPI_Barrier(MPI_COMM_WORLD); // Allows all other processes to start once user input has been received.

    string icPath = dir_path + "/Data/ic.txt";
    string finalFieldPath = dir_path + "/Data/finalField.txt";
    //string valsPerLoopPath = dir_path + "/Data/valsPerLoop_" + input + ".txt";
    // string testPath = dir_path + "/Data/test_" + to_string(rank) + ".txt";
    // string testMergePath = dir_path + "/Data/test_" + gifTag + ".txt";

    ifstream ic (icPath.c_str());
    ofstream finalField (finalFieldPath.c_str());
    //ofstream valsPerLoop (valsPerLoopPath.c_str());
    // ofstream test (testPath.c_str());
    // ofstream testMerge (testMergePath.c_str());

    x0 = int(0.5*(nx-1));
    y0 = int(0.5*(ny-1));
    z0 = int(0.5*(nz-1));

    int coreStart, coreEnd;
    if(rank < chunkRem){ coreStart = rank*(chunk+1); coreEnd = (rank+1)*(chunk+1); }
    else{ coreStart = rank*chunk+chunkRem; coreEnd = (rank+1)*chunk+chunkRem; }
    int dataStart = coreStart-haloSize;
    int dataEnd = coreEnd+haloSize;

    double wasteData[5];

    if(ic_type=="stationary data"){

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

    } else if(ic_type=="data"){

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

    } else if(ic_type=="random"){

    	// Use the seed to generate the data
		mt19937 generator (seed);
        normal_distribution<double> distribution (mean,stdev);

        // Skip the random numbers ahead to the appropriate point
        for(i=0;i<coreStart;i++){
        	for(comp=0;comp<2;comp++){ // Make sure to change this if randomly generating gauge fields too rather than leaving them at zero.

        		double randomWaste = distribution(generator);

        	}
        }



        for(i=haloSize;i<coreSize+haloSize;i++){

        	phi[i] = distribution(generator);
        	phi[totSize*nts+i] = distribution(generator);

        	// Set next timestep as equal to the first
        	phi[totSize+i] = phi[i];
        	phi[totSize*(nts+1)+i] = phi[totSize*nts+i];

        	// Leave the gauge fields as zero (set by initialisation)

        }

        //cout << "Rank " << rank << "has phi[haloSize] = " << phi[haloSize] << ", and the next random number would be " << distribution(generator) << endl;

        // Now that the core data has been generated, need to communicate the haloes between processes. 

        for(comp=0;comp<2;comp++){ 

        	MPI_Sendrecv(&phi[totSize*nts*comp+haloSize],haloSize,MPI_DOUBLE,(rank-1+size)%size,0, // Send this
        				 &phi[totSize*nts*comp+coreSize+haloSize],haloSize,MPI_DOUBLE,(rank+1)%size,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE); // Receive this

        	MPI_Sendrecv(&phi[totSize*nts*comp+coreSize],haloSize,MPI_DOUBLE,(rank+1)%size,0,
        				 &phi[totSize*nts*comp],haloSize,MPI_DOUBLE,(rank-1+size)%size,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

        	MPI_Sendrecv(&phi[totSize*(nts*comp+1)+haloSize],haloSize,MPI_DOUBLE,(rank-1+size)%size,1, // Send this
        				 &phi[totSize*(nts*comp+1)+coreSize+haloSize],haloSize,MPI_DOUBLE,(rank+1)%size,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE); // Receive this

        	MPI_Sendrecv(&phi[totSize*(nts*comp+1)+coreSize],haloSize,MPI_DOUBLE,(rank+1)%size,1,
        				 &phi[totSize*(nts*comp+1)],haloSize,MPI_DOUBLE,(rank-1+size)%size,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

        }

        // for(comp=0;comp<3;comp++){ 

        // 	MPI_Sendrecv(&theta[totSize*(nts*comp+tPast)+haloSize],haloSize,MPI_DOUBLE,(rank-1+size)%size,0,
        // 				 &theta[totSize*(nts*comp+tPast)+coreSize+haloSize],haloSize,MPI_DOUBLE,(rank+1)%size,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

        // 	MPI_Sendrecv(&theta[totSize*(nts*comp+tPast)+coreSize],haloSize,MPI_DOUBLE,(rank+1)%size,0,
        // 				 &theta[totSize*(nts*comp+tPast)],haloSize,MPI_DOUBLE,(rank-1+size)%size,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

        // }   	

    }

    gettimeofday(&end,NULL);

    if(rank==0){ cout << "Initial data loaded/generated in: " << end.tv_sec - start.tv_sec << "s" << endl; }

    gifStringPosFrame = 0;
    counter = 0;
    bool stringsExist = true;

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

        vector<double> xString, yString, zString; // Declares the vectors and clears them at every loop

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

        	// Calculate where the strings cross through grid faces. Each point searches the faces on the positive side (i.e pxpy, pxpz and pypz) to avoid double counting

	        if(stringDetect and (!detectBuffer or TimeStep>=damped_nt)){

	        	double x,y,z,coeff1[2],coeff2[2],coeff3[2],coeff4[2],a,b,c,discrim,sol1,sol2;
	        	int ipxpy,ipxpz,ipypz;

	        	// Need a few more indices and the physical coordinates in space

	        	ipxpy = ipy+ny*nz;
	        	ipxpz = ipz+ny*nz;
	        	ipypz = (ipy+dataStart+1)%nz + ( (ipy+dataStart)/nz )*nz - dataStart;

	        	x = ( (i+dataStart)/(ny*nz) - x0 )*dx;
	        	y = ( ((i+dataStart)/nz)%ny - y0 )*dy;
	        	z = ( (i+dataStart)%nz - z0 )*dz;


		        for(comp=0;comp<2;comp++){

		            // Do the same process for the real and imaginary components

		            coeff1[comp] = phi[totSize*(nts*comp+tNow)+ipxpy] - phi[totSize*(nts*comp+tNow)+ipx] - phi[totSize*(nts*comp+tNow)+ipy] + phi[totSize*(nts*comp+tNow)+i];
		            coeff2[comp] = (y+dy)*(phi[totSize*(nts*comp+tNow)+ipx] - phi[totSize*(nts*comp+tNow)+i]) + y*(phi[totSize*(nts*comp+tNow)+ipy] - phi[totSize*(nts*comp+tNow)+ipxpy]);
		            coeff3[comp] = (x+dx)*(phi[totSize*(nts*comp+tNow)+ipy] - phi[totSize*(nts*comp+tNow)+i]) + x*(phi[totSize*(nts*comp+tNow)+ipx] - phi[totSize*(nts*comp+tNow)+ipxpy]);
		            coeff4[comp] = (x+dx)*( (y+dy)*phi[totSize*(nts*comp+tNow)+i] - y*phi[totSize*(nts*comp+tNow)+ipy] ) - x*( (y+dy)*phi[totSize*(nts*comp+tNow)+ipx] - y*phi[totSize*(nts*comp+tNow)+ipxpy] );

		        }

		        // Substituting one equation into the other gives a quadratic equation for one of the coordinates. Now calculate the coefficients of the quadratic.

		        a = coeff1[1]*coeff3[0] - coeff1[0]*coeff3[1];  
		        b = coeff1[1]*coeff4[0] + coeff2[1]*coeff3[0] - coeff1[0]*coeff4[1] - coeff2[0]*coeff3[1];
		        c = coeff2[1]*coeff4[0] - coeff2[0]*coeff4[1];

		        discrim = b*b - 4*a*c;

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

		                        xString.push_back(sol1);
		                        yString.push_back(sol2);
		                        zString.push_back(z);

		                    }

		                }

		            }

		            // Otherwise no solutions

		        } else if(discrim >= 0){

		            // There will be two solutions (or a repeated, this may cause trouble). First solution is

		            sol2 = ( -b + sqrt(discrim) )/(2*a);

		            if(sol2>=y && sol2<=y+dy){

		                sol1 = -(coeff3[0]*sol2 + coeff4[0])/(coeff1[0]*sol2 + coeff2[0]);

		                if(sol1>=x && sol1<=x+dx){

		                    xString.push_back(sol1);
		                    yString.push_back(sol2);
		                    zString.push_back(z);

		                }

		            }

		            // Second solution is

		            sol2 = ( -b - sqrt(discrim) )/(2*a);

		            if(sol2>=y && sol2<=y+dy){

		                sol1 = -(coeff3[0]*sol2 + coeff4[0])/(coeff1[0]*sol2 + coeff2[0]);

		                if(sol1>=x && sol1<=x+dx){

		                    xString.push_back(sol1);
		                    yString.push_back(sol2);
		                    zString.push_back(z);

		                }

		            }

		        }

		        // Now repeat this process for the y directed face

	            for(comp=0;comp<2;comp++){

	                coeff1[comp] = phi[totSize*(nts*comp+tNow)+ipxpz] - phi[totSize*(nts*comp+tNow)+ipx] - phi[totSize*(nts*comp+tNow)+ipz] + phi[totSize*(nts*comp+tNow)+i];
	                coeff2[comp] = (z+dz)*(phi[totSize*(nts*comp+tNow)+ipx] - phi[totSize*(nts*comp+tNow)+i]) + z*(phi[totSize*(nts*comp+tNow)+ipz] - phi[totSize*(nts*comp+tNow)+ipxpz]);
	                coeff3[comp] = (x+dx)*(phi[totSize*(nts*comp+tNow)+ipz] - phi[totSize*(nts*comp+tNow)+i]) + x*(phi[totSize*(nts*comp+tNow)+ipx] - phi[totSize*(nts*comp+tNow)+ipxpz]);
	                coeff4[comp] = (x+dx)*( (z+dz)*phi[totSize*(nts*comp+tNow)+i] - z*phi[totSize*(nts*comp+tNow)+ipz] ) - x*( (z+dz)*phi[totSize*(nts*comp+tNow)+ipx] - z*phi[totSize*(nts*comp+tNow)+ipxpz] );

	            }

	            a = coeff1[1]*coeff3[0] - coeff1[0]*coeff3[1];
	            b = coeff1[1]*coeff4[0] + coeff2[1]*coeff3[0] - coeff1[0]*coeff4[1] - coeff2[0]*coeff3[1];
	            c = coeff2[1]*coeff4[0] - coeff2[0]*coeff4[1];

	            discrim = b*b - 4*a*c;

	            if(a==0){
	                if(b!=0){

	                    sol2 = -c/b;

	                    if(sol2>=z && sol2<=z+dz){

	                        sol1 = -(coeff3[0]*sol2 + coeff4[0])/(coeff1[0]*sol2 + coeff2[0]);

	                        if(sol1>=x && sol1<=x+dx){

	                            xString.push_back(sol1);
	                            yString.push_back(y);
	                            zString.push_back(sol2);

	                        }

	                    }

	                }
	            } else if(discrim >= 0){

	                sol2 = ( -b + sqrt(discrim) )/(2*a);

	                if(sol2>=z && sol2<=z+dz){

	                    sol1 = -(coeff3[0]*sol2 + coeff4[0])/(coeff1[0]*sol2 + coeff2[0]);

	                    if(sol1>=x && sol1<=x+dx){

	                        xString.push_back(sol1);
	                        yString.push_back(y);
	                        zString.push_back(sol2);

	                    }

	                }

	                sol2 = ( -b - sqrt(discrim) )/(2*a);

	                if(sol2>=z && sol2<=z+dz){

	                    sol1 = -(coeff3[0]*sol2 + coeff4[0])/(coeff1[0]*sol2 + coeff2[0]);

	                    if(sol1>=x && sol1<=x+dx){

	                        xString.push_back(sol1);
	                        yString.push_back(y);
	                        zString.push_back(sol2);

	                    }

	                }

	            }



	            // Now repeat one more time for the x directed face

	            for(comp=0;comp<2;comp++){

	                coeff1[comp] = phi[totSize*(nts*comp+tNow)+ipypz] - phi[totSize*(nts*comp+tNow)+ipy] - phi[totSize*(nts*comp+tNow)+ipz] + phi[totSize*(nts*comp+tNow)+i];
	                coeff2[comp] = (z+dz)*(phi[totSize*(nts*comp+tNow)+ipy] - phi[totSize*(nts*comp+tNow)+i]) + z*(phi[totSize*(nts*comp+tNow)+ipz] - phi[totSize*(nts*comp+tNow)+ipypz]);
	                coeff3[comp] = (y+dy)*(phi[totSize*(nts*comp+tNow)+ipz] - phi[totSize*(nts*comp+tNow)+i]) + y*(phi[totSize*(nts*comp+tNow)+ipy] - phi[totSize*(nts*comp+tNow)+ipypz]);
	                coeff4[comp] = (y+dy)*( (z+dz)*phi[totSize*(nts*comp+tNow)+i] - z*phi[totSize*(nts*comp+tNow)+ipz] ) - y*( (z+dz)*phi[totSize*(nts*comp+tNow)+ipy] - z*phi[totSize*(nts*comp+tNow)+ipypz] );

	            }

	            a = coeff1[1]*coeff3[0] - coeff1[0]*coeff3[1];
	            b = coeff1[1]*coeff4[0] + coeff2[1]*coeff3[0] - coeff1[0]*coeff4[1] - coeff2[0]*coeff3[1];
	            c = coeff2[1]*coeff4[0] - coeff2[0]*coeff4[1];

	            discrim = b*b - 4*a*c;

	            if(a==0){
	                if(b!=0){

	                    sol2 = -c/b;

	                    if(sol2>=z && sol2<=z+dz){

	                        sol1 = -(coeff3[0]*sol2 + coeff4[0])/(coeff1[0]*sol2 + coeff2[0]);

	                        if(sol1>=y && sol1<=y+dy){

	                            xString.push_back(x);
	                            yString.push_back(sol1);
	                            zString.push_back(sol2);

	                        }

	                    }

	                }
	            } else if(discrim >= 0){

	                sol2 = ( -b + sqrt(discrim) )/(2*a);

	                if(sol2>=z && sol2<=z+dz){

	                    sol1 = -(coeff3[0]*sol2 + coeff4[0])/(coeff1[0]*sol2 + coeff2[0]);

	                    if(sol1>=y && sol1<=y+dy){

	                        xString.push_back(x);
	                        yString.push_back(sol1);
	                        zString.push_back(sol2);

	                    }

	                }

	                sol2 = ( -b - sqrt(discrim) )/(2*a);

	                if(sol2>=z && sol2<=z+dz){

	                    sol1 = -(coeff3[0]*sol2 + coeff4[0])/(coeff1[0]*sol2 + coeff2[0]);

	                    if(sol1>=y && sol1<=y+dy){

	                        xString.push_back(x);
	                        yString.push_back(sol1);
	                        zString.push_back(sol2);

	                    }

	                }

	            }

		    }

        }

        // Next step is to collect all the intersections together on the master process and output them to text

        if(stringDetect and (!detectBuffer or TimeStep>=damped_nt)){

	        if(rank==0){

	        	vector<double> xStringFull, yStringFull, zStringFull;

	        	xStringFull.insert(xStringFull.end(),xString.begin(),xString.end());
	        	yStringFull.insert(yStringFull.end(),yString.begin(),yString.end());
	        	zStringFull.insert(zStringFull.end(),zString.begin(),zString.end());

	        	for(i=1;i<size;i++){

	        		// Find out how many intersections each process found

	        		int nIntersections;

	        		MPI_Recv(&nIntersections,1,MPI_INT,i,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

	        		// Declare vectors to receive each process's intersections and receive the data.

	        		vector<double> localxString(nIntersections,0.0), localyString(nIntersections,0.0), localzString(nIntersections,0.0);

	        		MPI_Recv(&localxString[0],nIntersections,MPI_DOUBLE,i,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	        		MPI_Recv(&localyString[0],nIntersections,MPI_DOUBLE,i,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	        		MPI_Recv(&localzString[0],nIntersections,MPI_DOUBLE,i,4,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

	        		// Append the data to the full collection of intersections

	        		xStringFull.insert(xStringFull.end(),localxString.begin(),localxString.end());
	        		yStringFull.insert(yStringFull.end(),localyString.begin(),localyString.end());
	        		zStringFull.insert(zStringFull.end(),localzString.begin(),localzString.end());

	        	}

	        	// Full set of intersections are now collected on the master process. Now just output them to text.

	        	if(makeStringPosGif and TimeStep%saveFreq == 0 and stringsExist){

	                ss.str(string());
	                ss << gifStringPosFrame;
	                string gifStringPosDataPath = dir_path + "/GifData/gifStringPosData_" + gifTag + "_nx_" + to_string(nx) + "_seed_" + to_string(round(seed)) + ss.str() + ".txt";
	                ofstream gifStringPosData (gifStringPosDataPath.c_str());
	                gifStringPosFrame+=1;

	                for(i=0;i<xStringFull.size();i++){

	                    gifStringPosData << xStringFull[i] << " " << yStringFull[i] << " " << zStringFull[i] << endl;

	                }

	                if(xStringFull.size() == 0){ stringsExist = false; }
	            }

	        } else{

	        	// Tell the master process how many string intersections to expect
	        	int nIntersections = xString.size();
	        	MPI_Send(&nIntersections,1,MPI_INT,0,1,MPI_COMM_WORLD);

	        	// Send the data
	        	MPI_Send(&xString[0],nIntersections,MPI_DOUBLE,0,2,MPI_COMM_WORLD);
	        	MPI_Send(&yString[0],nIntersections,MPI_DOUBLE,0,3,MPI_COMM_WORLD);
	        	MPI_Send(&zString[0],nIntersections,MPI_DOUBLE,0,4,MPI_COMM_WORLD);

	        }

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

        	MPI_Sendrecv(&phi[totSize*(nts*comp+tPast)+haloSize],haloSize,MPI_DOUBLE,(rank-1+size)%size,0, // Send this
        				 &phi[totSize*(nts*comp+tPast)+coreSize+haloSize],haloSize,MPI_DOUBLE,(rank+1)%size,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE); // Receive this

        	MPI_Sendrecv(&phi[totSize*(nts*comp+tPast)+coreSize],haloSize,MPI_DOUBLE,(rank+1)%size,0,
        				 &phi[totSize*(nts*comp+tPast)],haloSize,MPI_DOUBLE,(rank-1+size)%size,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

        }

        for(comp=0;comp<3;comp++){ 

        	MPI_Sendrecv(&theta[totSize*(nts*comp+tPast)+haloSize],haloSize,MPI_DOUBLE,(rank-1+size)%size,0,
        				 &theta[totSize*(nts*comp+tPast)+coreSize+haloSize],haloSize,MPI_DOUBLE,(rank+1)%size,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

        	MPI_Sendrecv(&theta[totSize*(nts*comp+tPast)+coreSize],haloSize,MPI_DOUBLE,(rank+1)%size,0,
        				 &theta[totSize*(nts*comp+tPast)],haloSize,MPI_DOUBLE,(rank-1+size)%size,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

        }

        // Code for testing the output

  //       if(TimeStep==nt-1){

	 //        if(rank==0){

	 //        	vector<double> phittOut(2*nPos,0.0), thetattOut(3*nPos,0.0);
	 //        	for(i=0;i<coreSize;i++){ phittOut[i] = phitt[i]; phittOut[nPos+i] = phitt[coreSize+i]; thetattOut[i] = thetatt[i]; thetattOut[nPos+i] = thetatt[coreSize+i]; thetattOut[2*nPos+i] = thetatt[2*coreSize+i]; }

	 //        	for(i=1;i<size;i++){

	 //        		int localCoreStart;
	 //        		int localCoreSize;
	 //        		if(i<chunkRem){ localCoreStart = i*(chunk+1); localCoreSize = chunk+1; }
	 //        		else{ localCoreStart = i*chunk + chunkRem; localCoreSize = chunk; }

	 //        		MPI_Recv(&phittOut[localCoreStart],localCoreSize,MPI_DOUBLE,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	 //        		MPI_Recv(&phittOut[nPos+localCoreStart],localCoreSize,MPI_DOUBLE,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

	 //        		MPI_Recv(&thetattOut[localCoreStart],localCoreSize,MPI_DOUBLE,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	 //        		MPI_Recv(&thetattOut[nPos+localCoreStart],localCoreSize,MPI_DOUBLE,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	 //        		MPI_Recv(&thetattOut[2*nPos+localCoreStart],localCoreSize,MPI_DOUBLE,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);  

		// 		}

	 //        	//for(i=0;i<nPos;i++){ testMerge << phittOut[i] << " " << phittOut[nPos+i] << " " << thetattOut[i] << " " << thetattOut[nPos+i] << " " << thetattOut[2*nPos+i] << endl; }

	 //        	// testArray[1] = 10;

	 //        	// cout << *foo << " " << *foo2 << endl;


	 //        } else{ MPI_Send(&phitt[0],coreSize,MPI_DOUBLE,0,0,MPI_COMM_WORLD); MPI_Send(&phitt[coreSize],coreSize,MPI_DOUBLE,0,0,MPI_COMM_WORLD); MPI_Send(&thetatt[0],coreSize,MPI_DOUBLE,0,0,MPI_COMM_WORLD); MPI_Send(&thetatt[coreSize],coreSize,MPI_DOUBLE,0,0,MPI_COMM_WORLD); MPI_Send(&thetatt[2*coreSize],coreSize,MPI_DOUBLE,0,0,MPI_COMM_WORLD); }

		// }

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