//intialiseSurface.cpp: initialise the surface geometries
//Marco, 29 June 2011

#include "wet.h"
#include <cmath>
using namespace std;

void wet::initialiseSurface(void)
{
	computeCoordinate(k);
	if( dimensions==2)
    {
        if (xk>Px and xk<Px+Dx and yk>Py and yk<Py+Dy)
        {mask[k]=28;}
    }


    if( dimensions==3)
    {
        if (xk>Px and xk<Px+Dx and yk>Py and yk<Py+Dy and zk>Pz and zk<Pz+Dh)
        {mask[k]=28;}
        
        
    }





}

