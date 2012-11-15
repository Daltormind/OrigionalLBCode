//LBAlgorithm.cpp: lattice boltzmann algorithm
//Marco, 30 June 2011

#include "wet.h"

void wet::propagation(void)
{
	for( k = k1 ; k < k2 ; k++)
	{
		if (mask[k] != 28) 
		{
			ff0[k] = fn0[k]; 
			ff1[dd1[k]] = fn1[k]; ff2[dd2[k]] = fn2[k]; ff3[dd3[k]] = fn3[k]; ff4[dd4[k]] = fn4[k]; ff5[dd5[k]] = fn5[k]; ff6[dd6[k]] = fn6[k];
			ffa[dda[k]] = fna[k]; ffb[ddb[k]] = fnb[k]; ffc[ddc[k]] = fnc[k]; ffd[ddd[k]] = fnd[k]; ffe[dde[k]] = fne[k]; fff[ddf[k]] = fnf[k]; 
			ffg[ddg[k]] = fng[k]; ffh[ddh[k]] = fnh[k]; ffi[ddi[k]] = fni[k]; ffj[ddj[k]] = fnj[k]; ffk[ddk[k]] = fnk[k]; ffl[ddl[k]] = fnl[k];

			gg0[k] = gn0[k]; 
			gg1[dd1[k]] = gn1[k]; gg2[dd2[k]] = gn2[k]; gg3[dd3[k]] = gn3[k]; gg4[dd4[k]] = gn4[k]; gg5[dd5[k]] = gn5[k]; gg6[dd6[k]] = gn6[k];
			gga[dda[k]] = gna[k]; ggb[ddb[k]] = gnb[k]; ggc[ddc[k]] = gnc[k]; ggd[ddd[k]] = gnd[k]; gge[dde[k]] = gne[k]; ggf[ddf[k]] = gnf[k];
			ggg[ddg[k]] = gng[k]; ggh[ddh[k]] = gnh[k]; ggi[ddi[k]] = gni[k]; ggj[ddj[k]] = gnj[k]; ggk[ddk[k]] = gnk[k]; ggl[ddl[k]] = gnl[k];
		}
	}
}
