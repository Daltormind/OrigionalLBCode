//ComputeCoordinate.cpp: This function computes the 3-D coordinate of index k. Filling order: zk, yk, xk 
//Marco, 29 June 2011

#include "wet.h"

using namespace std;

void wet::computeCoordinate(int i)
{  
	int kTot;

	if(rank<ROOT)
	{
		kTot=i+rank*(k2-k1)-k1;
		xk=int(kTot/(float) (LZ*LY));		//Index node of x coordinate 
		yk=int((kTot-xk*LZ*LY)/(float) LZ);	//Index node of y coordinate
		zk=kTot-xk*LZ*LY-yk*LZ;			//Index node of z coordinate
	}
	if(rank==ROOT)
	{
		kTot=i+rank*(k2-(LX%size+1)*k1)-k1;
		xk=int(kTot/(float) (LZ*LY));		//Index node of x coordinate 
		yk=int((kTot-xk*LZ*LY)/(float) LZ);	//Index node of y coordinate
		zk=kTot-xk*LZ*LY-yk*LZ;			//Index node of z coordinate
	}
	if(rank>ROOT)
	{
		kTot=i+rank*(k2-k1)+(LX%size-1)*k1;
		xk=int(kTot/(float) (LZ*LY));		//Index node of x coordinate 
		yk=int((kTot-xk*LZ*LY)/(float) LZ);	//Index node of y coordinate
		zk=kTot-xk*LZ*LY-yk*LZ;			//Index node of z coordinate
	}

}
