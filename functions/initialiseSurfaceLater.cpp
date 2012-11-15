//intialiseSurfaceLater.cpp: initialise the surface geometries, and includes setting up mask
//necessary if the surface should be included after drop equilibration for teta_CB_fix=180 deg
//Lisa

#include "wet.h"
#include <cmath>
using namespace std;

void wet::initialiseSurfaceLater(void)
{	
	//cout << "** ** Process " << rank << ": initSurfaceLater started " << endl;
	for (k = k1; k < k2; k++) 
		
		initialiseSurface(); //initialise surface drawing mask=28 if the point xk, yk, zk is part of solid surface
	//cout << "** ** Process " << rank << ": initSurfaceLater check1 " << endl;
	exchangeMask();
	//cout << "** ** Process " << rank << ": initSurfaceLater check2 " << endl;
	for (k = k1; k < k2; k++)
		relabel();		//Redrawing mask to put value 1 near the surface
	
	for (k = k1; k < k2; k++) 
		leakSearch();		//Search for missing point in the relabel function 
	
	//cout << "** ** Process " << rank << ": initSurfaceLater check3 " << endl;
	generateGlobalMask();
	
	//cout << "** ** Process " << rank << ": initSurfaceLater check4 " << endl;
	
	for (k = k1; k < k2; k++) {
		
		
		if (mask[k] == 28) 
			
		{	//cout << "-- Process " << rank << " k = " << k << " mask = " << mask[k] << endl;
			n[k] = 1.0; 
			p[k] = 0.0;
			uxs[k] = 0.0;
			uys[k] = 0.0;
			uzs[k] = 0.0;
		}
	}
	//cout << "** ** Process " << rank << ": initSurfaceLater check5 " << endl;
}	

