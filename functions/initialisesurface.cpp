//intialiseSurface.cpp: initialise the surface geometries
//Marco, 29 June 2011

#include "wet.h"
#include <cmath>
using namespace std;

void wet::initialiseSurface(void)
{
	computeCoordinate(k);
	
	if (xk>Px and xk<Px+Dx and yk>Py and yk<Py+Dy)
    {mask[k]=28;}
}	

