#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <sys/time.h>
#include <omp.h>
#include <complex>

#include "header.hpp"

using namespace std;

// Constructs the initial field configuration for a completely straight string by interpolating from radial profiles.
// Used to test the evolution code 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                 		  Parameters & Declarations
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

const int nx = 101;
const int ny = 101;
const int nz = 51;
const double dx = 0.5;
const double dy = 0.5;
const double dz = 0.5;

const int n = 1;
const double g = 0.5;

const int SORnx = 20001;
const double SORa = 0.01;

int main(){

	Array phi(2,nx,ny,nz,0.0), A(3,nx,ny,nz,0.0), SOR_Fields(SORnx,2,0.0);

	double x,y,z,distance,phiMag,AMag;

	int i,j,k,pClosest;

	struct timeval start, end;
    gettimeofday(&start, NULL);

    string file_path = __FILE__;
    string dir_path = file_path.substr(0,file_path.find_last_of('/'));

    string SOR_inputPath = dir_path + "/SOR_Fields.txt";
    string icPath = dir_path + "/ic.txt";

    ifstream SOR_input (SOR_inputPath.c_str());
    ofstream ic (icPath.c_str());

    for(i=0;i<SORnx;i++){

    	SOR_input >> SOR_Fields(i,0) >> SOR_Fields(i,1);

    }

    int x0 = round((nx-1)/2);
    int y0 = round((ny-1)/2);
    int z0 = round((nz-1)/2);

	for(i=0;i<nx;i++){

		x = (i-x0)*dx;

		for(j=0;j<ny;j++){

			y = (j-y0)*dy;

			for(k=0;k<nz;k++){

				distance = sqrt(x*x + y*y);

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

     			if(x==0){

     				// To prevent division by zero

     				phi(0,i,j,k) = 0;
     				A(1,i,j,k) = 0;

     			} else{

     				phi(0,i,j,k) = phiMag*x/distance;
     				A(1,i,j,k) = AMag*x/pow(distance,2);

     			}

     			if(y==0){

     				// To prevent division by zero

     				phi(1,i,j,k) = 0;
     				A(0,i,j,k) = 0;

     			} else{

     				phi(1,i,j,k) = phiMag*y/distance;
     				A(0,i,j,k) = -AMag*y/(pow(distance,2));

     			}

     			A(2,i,j,k) = 0;

     		}
		}
	}

	for(i=0;i<nx;i++){
    	for(j=0;j<ny;j++){
    		for(k=0;k<nz;k++){

    			// Convert gauge field to lattice link variable for use in evolution code

    			ic << phi(0,i,j,k) << " " << phi(1,i,j,k) << " " << dx*g*A(0,i,j,k) << " " << dy*g*A(1,i,j,k) << " " << dz*g*A(2,i,j,k) << endl;

    		}
    	}
    }

}