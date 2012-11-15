// computeContactArea.cpp: calculates the area of contact between water and substrate, 
// as well as the the length (actually the area) of the contact line with two definitions
// Lisa, 12 Aug 2011

#include "wet.h"
#include <cmath>
using namespace std;

void wet::computeContactArea(void)
{
	surfArea = 0;
	clLength = 0;
	clLength_n = 0;
	
	LFree = 0;

	for (k = k1; k < k2 ; k++){
		if (mask[k] == 1) {
			if (p[k] >= 0.0) surfArea += p[k]; 
			
			computeCoordinate(k);
			if (p[k] < -tanh(1.0) && yk == 0) LFree += 2; //count nodes that are: 1. gas (not interface) 2. edge (y=0 and surface layer); x2 for both sides of ridge
			//cout << "Process " << rank << ": , t = " << t << " in computeContactArea, k = " << k << " x " << xk << " y " << yk << " z " << zk << ": LFree = "<< LFree << endl;
			
			if (fabs(p[k]) <= p_thresh) {
				clLength += 1-fabs(p[k]);  
				
			}
			if (fabs(p[k]) <= tanh(1.0)) { 
				clLength_n += 1-fabs(p[k]);
				//clLength_n2 += 1;
			}
		}
		
	}
	double reducedValue;
	
	reducedValue = 0.0;
	MPI_Reduce(&surfArea,&reducedValue,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	surfArea = reducedValue;
	
	reducedValue = 0.0;
	MPI_Reduce(&LFree,&reducedValue,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	LFree = reducedValue;
	
	
	reducedValue = 0.0;
	MPI_Reduce(&clLength,&reducedValue,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	clLength = reducedValue;
	
	reducedValue = 0.0;
	MPI_Reduce(&clLength_n,&reducedValue,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	clLength_n = reducedValue;
	
}






