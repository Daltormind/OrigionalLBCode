//computeFreeEnergy.cpp: calculates t
//Marco, 30 June 2011

#include "wet.h"
#include <cmath>
using namespace std;

void wet::computeFreeEnergy(void)   //corrected calculation, called in equilibrium
{   
    
    //if (mask[k] != 1) {   // only for testing exclusion of some areas in dissipation calculation (ARolling, computeDissipation, computeFreeEnergy)
//computeCoordinate(k);
//if ( zk > (Dh+dropletCenterZ) && (zk + dropletCenterZ) < LZ){	
		
	double gradsq = dp_dx*dp_dx + dp_dy*dp_dy + dp_dz*dp_dz;
    
	double en = A*(-0.5*p[k]*p[k] + 0.25*p[k]*p[k]*p[k]*p[k]) + kappa_p/2.0*gradsq - A*(-0.5 + 0.25);
    
	//for set interface<->bulk threshold p_thresh:
	//
	//if (rank == ROOT) cout << "* * * Process " << rank << ", k = " << k << " compFreeE check1" << endl;
	
	if (fabs(p[k]) > p_thresh)
		bulkE += c2/3.0*n[k]*log(n[k]) + en;
	
	if (fabs(p[k]) <= p_thresh)
		bulkE += c2/3.0*n[k]*log(n[k]);

	if (fabs(p[k]) <= p_thresh && mask[k] != 1)
		interfaceE += en;
			
	if (mask[k] == 1){
		surfaceE += phi11*kappa_p*p[k];    
		if (fabs(p[k]) <= p_thresh) 
			surfaceE += en;
	}
		
	//for "natural threshold" p_thresh_n = tanh(1):		
	
	if (fabs(p[k]) > p_thresh_n){
		bulkE_n += c2/3.0*n[k]*log(n[k]) + en;
		//eBal[k]+= c2/3.0*n[k]*log(n[k]) + en;
	}
	if (fabs(p[k]) <= p_thresh_n){
		bulkE_n += c2/3.0*n[k]*log(n[k]);
		//eBal[k] += c2/3.0*n[k]*log(n[k]);
	}
	if (mask[k] != 1 && fabs(p[k]) <= p_thresh_n){
		interfaceE_n += en;
		//eBal[k]      += en;
	}
	if (mask[k] == 1){
		surfaceE_n += phi11*kappa_p*p[k];  
		//eBal[k]    += phi11*kappa_p*p[k];
		if (fabs(p[k]) <= p_thresh_n){
			surfaceE_n += en; 
			//eBal[k]    += en;
		}
	}
	//energy_n = bulkE_n + interfaceE_n + surfaceE_n;
	
		
	//}
	
	if (k==(k2-1)){
		double reducedEnergy;
		
		if(t%(infoStep) == 0){
			
			reducedEnergy = 0.0;
			MPI_Reduce(&bulkE,&reducedEnergy,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
			bulkE = reducedEnergy;
						
			reducedEnergy = 0.0;
			MPI_Reduce(&interfaceE,&reducedEnergy,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
			interfaceE = reducedEnergy;
			
			reducedEnergy = 0.0;
			MPI_Reduce(&surfaceE,&reducedEnergy,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
			surfaceE = reducedEnergy;
			
			//energy = bulkE + interfaceE + surfaceE;
			
			reducedEnergy = 0.0;
			MPI_Reduce(&bulkE_n,&reducedEnergy,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
			bulkE_n = reducedEnergy;
			
			reducedEnergy = 0.0;
			MPI_Reduce(&interfaceE_n,&reducedEnergy,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
			interfaceE_n = reducedEnergy;
			
			reducedEnergy = 0.0;
			MPI_Reduce(&surfaceE_n,&reducedEnergy,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
			surfaceE_n = reducedEnergy;
			
			//energy_n = bulkE_n + interfaceE_n + surfaceE_n;
		}
		
		if (rank == ROOT) {
			
			energy = bulkE + interfaceE + surfaceE;
			energy_n = bulkE_n + interfaceE_n + surfaceE_n;
		}
	}
}




