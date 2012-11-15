//ComputeMomenta.cpp: this function computes n (density), phi (energy) and the velocity
//Marco, 29 June 2011

#include "wet.h"

void wet::computeMomenta(void)
{
	//if(t%infoStep==0)
	//	cout << "Process " << rank << ": computing momenta..." << endl;
	for( k = k1 ; k < k2 ; k++)
	{         
		if(mask[k]!=28)
		{                                                         // Works out the moments
		    
			n[k] = ff0[k] + ff1[k] + ff2[k] + ff3[k] + ff4[k] + ff5[k] + ff6[k] + ffa[k] + ffb[k] + ffc[k] + ffd[k] + ffe[k] + fff[k] + ffg[k] + ffh[k] + ffi[k] + ffj[k] + ffk[k] + ffl[k];
			
			uxs[k] = (ff1[k] + ffa[k] + ffc[k] + ffi[k] +ffk[k] - ff2[k] - ffb[k] - ffd[k] - ffj[k] - ffl[k])/n[k];
			uys[k] = (ff3[k] + ffa[k] + ffb[k] + ffe[k] + ffg[k] - ff4[k] - ffc[k] - ffd[k] - fff[k] - ffh[k])/n[k];
			uzs[k] = (ff5[k] + ffe[k] + fff[k] + ffi[k] + ffj[k] - ff6[k] - ffg[k] - ffh[k] - ffk[k] - ffl[k])/n[k];

			p[k] = gg0[k]+gg1[k]+gg2[k]+gg3[k]+gg4[k]+gg5[k]+gg6[k]+gga[k]+ggb[k]+ggc[k]+ggd[k]+gge[k]+ggf[k]+
			ggg[k]+ggh[k]+ggi[k]+ggj[k]+ggk[k]+ggl[k];

		
		}
	}

	//if(t%infoStep==0)
	//	cout << "Process " << rank << ": momenta computed." << endl;	
}
