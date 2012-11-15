//LBAlgorithm.cpp: lattice boltzmann algorithm
//Marco, 30 June 2011

#include "wet.h"

void wet::propagation(void)
{
	
	for( k = k1; k < k2; k++)
	{
		if (mask[k] != 28) 
		{

			//Propagation of densities non moving along x
			ff0[k] = fn0[k]; 
			ff3[k] = fn3[dd4[k]];
			ff4[k] = fn4[dd3[k]];
			ff5[k] = fn5[dd6[k]];
			ff6[k] = fn6[dd5[k]];
			ffe[k] = fne[ddh[k]]; 
			fff[k] = fnf[ddg[k]]; 
			ffg[k] = fng[ddf[k]]; 
			ffh[k] = fnh[dde[k]];


			gg0[k] = gn0[k];
			gg3[k] = gn3[dd4[k]]; 
			gg4[k] = gn4[dd3[k]]; 
			gg5[k] = gn5[dd6[k]]; 
			gg6[k] = gn6[dd5[k]];
			gge[k] = gne[ddh[k]];
			ggf[k] = gnf[ddg[k]];
			ggg[k] = gng[ddf[k]]; 
			ggh[k] = gnh[dde[k]];


			//Propagation of densities moving along x
			ff1[k] = fn1[dd2[k]];
			ffa[k] = fna[ddd[k]];
			ffc[k] = fnc[ddb[k]];
			ffi[k] = fni[ddl[k]];
			ffk[k] = fnk[ddj[k]];

			gg1[k] = gn1[dd2[k]];
			gga[k] = gna[ddd[k]];
			ggc[k] = gnc[ddb[k]]; 
			ggi[k] = gni[ddl[k]];
			ggk[k] = gnk[ddj[k]];

			//Propagation of densities moving along -x
 			ff2[k] = fn2[dd1[k]];   
			ffb[k] = fnb[ddc[k]];  
			ffd[k] = fnd[dda[k]];   
			ffj[k] = fnj[ddk[k]];  
			ffl[k] = fnl[ddi[k]];

			 
			gg2[k] = gn2[dd1[k]]; 
			ggb[k] = gnb[ddc[k]];  
			ggd[k] = gnd[dda[k]];   
			ggj[k] = gnj[ddk[k]];  
			ggl[k] = gnl[ddi[k]];
		}
	}
}
