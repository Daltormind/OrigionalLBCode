//intialiseSurface.cpp: initialise the surface geometries
//Marco, 29 June 2011

#include "wet.h"
#include <cmath>
using namespace std;

void wet::initialiseSurface(void)
{
	computeCoordinate(k);
	double xt;
	double yt;
	
	xt=abs(xk-Px);
	yt=abs(yk-Py);
	
	
	if( dimensions==2)
    {
        if ((xk>Px-Dx and xk<Px+Dx) or (xk>Py-Dx and xk<Py+Dx) )
        {mask[k]=28;}
    }


    if( dimensions==3)
    {
        if (xt<Dx and yt<Dy)
        {mask[k]=28;}
        
        
    }

	xt=abs(xk-PeriodX);
	yt=abs(yk-PeriodY);
	
	if( dimensions==3)
    {
        if (xt<Dx and yt<Dy)
        {mask[k]=28;}
        
        
    }



}

