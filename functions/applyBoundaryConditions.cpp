//LBAlgorithm.cpp: lattice boltzmann algorithm
//Marco, 30 June 2011

#include "wet.h"

void wet::applyBoundaryConditions(void)
{
  
	for(k = k1 ; k < k2 ; k++)
	{
		switch(mask[k]) 
		{ 
			case 0:
				break;
			case 1:
								// Bortolo's version:
				/*
				if (mask[dd1[k]] == 28)
				 {       //solid wall is on th mask[k]=1 lattice nodes
				 ff2[k] = ff1[k]; gg2[k] = gg1[k];
				 ff0[k] += fn1[k] - ff2[k]; gg0[k] += gn1[k] - gg2[k];   
				 //ff0[k] = fn1[k] - ff2[k]; gg0[k] = gn1[k] - gg2[k];
				 
				 }
				 if (mask[dd2[k]] == 28)
				 {
				 ff1[k] = ff2[k]; gg1[k] = gg2[k];
				 ff0[k] += fn2[k] - ff1[k]; gg0[k] += gn2[k] - gg1[k];  
				 //ff0[k] = fn2[k] - ff1[k]; gg0[k] = gn2[k] - gg1[k];
				 }
				 if (mask[dd3[k]] == 28)
				 {
				 ff4[k] = ff3[k];gg4[k] = gg3[k];
				 ff0[k] += fn3[k] - ff4[k]; gg0[k] += gn3[k] - gg4[k];
				 //ff0[k] = fn3[k] - ff4[k]; gg0[k] = gn3[k] - gg4[k];
				 }
				 if (mask[dd4[k]] == 28)
				 {
				 ff3[k] = ff4[k]; gg3[k] = gg4[k];
				 ff0[k] += fn4[k] - ff3[k]; gg0[k] += gn4[k] - gg3[k];
				 //ff0[k] = fn4[k] - ff3[k]; gg0[k] = gn4[k] - gg3[k];
				 }
				 if (mask[dd5[k]] == 28)
				 {
				 ff6[k] = ff5[k]; gg6[k] = gg5[k];
				 ff0[k] += fn5[k] - ff6[k]; gg0[k] += gn5[k] - gg6[k];
				 //ff0[k] = fn5[k] - ff6[k]; gg0[k] = gn5[k] - gg6[k];
				 }
				 if (mask[dd6[k]] == 28)
				 {
				 ff5[k] = ff6[k]; gg5[k] = gg6[k];
				 ff0[k] += fn6[k] - ff5[k]; gg0[k] += gn6[k] - gg5[k];
				 //ff0[k] = fn6[k] - ff5[k]; gg0[k] = gn6[k] - gg5[k];
				 }
				 
				 if (mask[dda[k]] == 28)
				 {
				 ffd[k] = ffa[k]; ggd[k] = gga[k];
				 ff0[k] += fna[k] - ffd[k]; gg0[k] += gna[k] - ggd[k]; 
				 //ff0[k] = fna[k] - ffd[k]; gg0[k] = gna[k] - ggd[k];
				 }
				 if (mask[ddb[k]] == 28)
				 {
				 ffc[k] = ffb[k]; ggc[k] = ggb[k];
				 ff0[k] += fnb[k] - ffc[k]; gg0[k] += gnb[k] - ggc[k];
				 //ff0[k] = fnb[k] - ffc[k]; gg0[k] = gnb[k] - ggc[k];
				 }
				 if (mask[ddc[k]] == 28)
				 {
				 ffb[k] = ffc[k]; ggb[k] = ggc[k];
				 ff0[k] += fnc[k] - ffb[k]; gg0[k] += gnc[k] - ggb[k]; 
				 //ff0[k] = fnc[k] - ffb[k]; gg0[k] = gnc[k] - ggb[k];
				 }
				 if (mask[ddd[k]] == 28)
				 {	
				 ffa[k] = ffd[k]; gga[k] = ggd[k];
				 ff0[k] += fnd[k] - ffa[k]; gg0[k] += gnd[k] - gga[k]; 
				 //ff0[k] = fnd[k] - ffa[k]; gg0[k] = gnd[k] - gga[k];
				 }
				 if (mask[dde[k]] == 28)
				 {
				 ffh[k] = ffe[k]; ggh[k] = gge[k];
				 ff0[k] += fne[k] - ffh[k]; gg0[k] += gne[k] - ggh[k];   
				 //ff0[k] = fne[k] - ffh[k]; gg0[k] = gne[k] - ggh[k];
				 }
				 
				 if (mask[ddf[k]] == 28)
				 {
				 ffg[k] = fff[k]; ggg[k] = ggf[k];
				 ff0[k] += fnf[k] - ffg[k]; gg0[k] += gnf[k] - ggg[k]; 
				 //ff0[k] = fnf[k] - ffg[k]; gg0[k] = gnf[k] - ggg[k]; 
				 }
				 if (mask[ddg[k]] == 28)
				 {
				 fff[k] = ffg[k]; ggf[k] = ggg[k];
				 ff0[k] += fng[k] - fff[k]; gg0[k] += gng[k] - ggf[k];
				 //ff0[k] = fng[k] - fff[k]; gg0[k] = gng[k] - ggf[k];
				 
				 }
				 if (mask[ddh[k]] == 28)
				 {
				 ffe[k] = ffh[k]; gge[k] = ggh[k];
				 ff0[k] += fnh[k] - ffe[k]; gg0[k] += gnh[k] - gge[k]; 
				 //ff0[k] = fnh[k] - ffe[k]; gg0[k] = gnh[k] - gge[k];  
				 }
				 if (mask[ddi[k]] == 28)
				 {
				 ffl[k] = ffi[k]; ggl[k] = ggi[k];
				 ff0[k] += fni[k] - ffl[k]; gg0[k] += gni[k] - ggl[k]; 
				 //ff0[k] = fni[k] - ffl[k]; gg0[k] = gni[k] - ggl[k];
				 }
				 if (mask[ddj[k]] == 28)
				 {
				 ffk[k] = ffj[k]; ggk[k] = ggj[k];
				 ff0[k] += fnj[k] - ffk[k]; gg0[k] += gnj[k] - ggk[k]; 
				 //ff0[k] = fnj[k] - ffk[k]; gg0[k] = gnj[k] - ggk[k];
				 }
				 if (mask[ddk[k]] == 28)
				 {
				 ffj[k] = ffk[k]; ggj[k] = ggk[k];
				 ff0[k] += fnk[k] - ffj[k]; gg0[k] += gnk[k] - ggj[k];  
				 //ff0[k] = fnk[k] - ffj[k]; gg0[k] = gnk[k] - ggj[k];
				 }	
				 if (mask[ddl[k]] == 28)
				 {
				 ffi[k] = ffl[k]; ggi[k] = ggl[k];
				 ff0[k] += fnl[k] - ffi[k]; gg0[k] += gnl[k] - ggi[k];  
				 //ff0[k] = fnl[k] - ffi[k]; gg0[k] = gnl[k] - ggi[k]; 
				 }
				*/ 
				 
				
				// Bortolo's version, Lisa-adjusted:
				if (mask[dd1[k]] == 28 || mask[dd1[k]] == 1)
				{       //solid wall is on th mask[k]=1 lattice nodes
					ff2[k] = ff1[k]; gg2[k] = gg1[k];
					ff0[k] += fn1[k] - ff2[k]; gg0[k] += gn1[k] - gg2[k];   
					//ff0[k] = fn1[k] - ff2[k]; gg0[k] = gn1[k] - gg2[k];

								}
				if (mask[dd2[k]] == 28 || mask[dd2[k]] == 1)
				{
					ff1[k] = ff2[k]; gg1[k] = gg2[k];
					ff0[k] += fn2[k] - ff1[k]; gg0[k] += gn2[k] - gg1[k];  
					//ff0[k] = fn2[k] - ff1[k]; gg0[k] = gn2[k] - gg1[k];
				}
				if (mask[dd3[k]] == 28 || mask[dd3[k]] == 1)
				{
					ff4[k] = ff3[k];gg4[k] = gg3[k];
					ff0[k] += fn3[k] - ff4[k]; gg0[k] += gn3[k] - gg4[k];
					//ff0[k] = fn3[k] - ff4[k]; gg0[k] = gn3[k] - gg4[k];
				}
				if (mask[dd4[k]] == 28 || mask[dd4[k]] == 1)
				{
					ff3[k] = ff4[k]; gg3[k] = gg4[k];
					ff0[k] += fn4[k] - ff3[k]; gg0[k] += gn4[k] - gg3[k];
					//ff0[k] = fn4[k] - ff3[k]; gg0[k] = gn4[k] - gg3[k];
				}
				if (mask[dd5[k]] == 28 || mask[dd5[k]] == 1)
				{
					ff6[k] = ff5[k]; gg6[k] = gg5[k];
					ff0[k] += fn5[k] - ff6[k]; gg0[k] += gn5[k] - gg6[k];
					//ff0[k] = fn5[k] - ff6[k]; gg0[k] = gn5[k] - gg6[k];
				}
				if ( mask[dd6[k]] == 28 || mask[dd6[k]] == 1)
				{
					ff5[k] = ff6[k]; gg5[k] = gg6[k];
					ff0[k] += fn6[k] - ff5[k]; gg0[k] += gn6[k] - gg5[k];
					//ff0[k] = fn6[k] - ff5[k]; gg0[k] = gn6[k] - gg5[k];
				}

				if (mask[dda[k]] == 28 || mask[dda[k]] == 1)
				{
					ffd[k] = ffa[k]; ggd[k] = gga[k];
					ff0[k] += fna[k] - ffd[k]; gg0[k] += gna[k] - ggd[k]; 
					//ff0[k] = fna[k] - ffd[k]; gg0[k] = gna[k] - ggd[k];
				}
				if (mask[ddb[k]] == 28 || mask[ddb[k]] == 1)
				{
					ffc[k] = ffb[k]; ggc[k] = ggb[k];
					ff0[k] += fnb[k] - ffc[k]; gg0[k] += gnb[k] - ggc[k];
					//ff0[k] = fnb[k] - ffc[k]; gg0[k] = gnb[k] - ggc[k];
				}
				if (mask[ddc[k]] == 28 || mask[ddc[k]] == 1)
				{
					ffb[k] = ffc[k]; ggb[k] = ggc[k];
					ff0[k] += fnc[k] - ffb[k]; gg0[k] += gnc[k] - ggb[k]; 
					//ff0[k] = fnc[k] - ffb[k]; gg0[k] = gnc[k] - ggb[k];
				}
				if (mask[ddd[k]] == 28 || mask[ddd[k]] == 1)
				{	
					ffa[k] = ffd[k]; gga[k] = ggd[k];
					ff0[k] += fnd[k] - ffa[k]; gg0[k] += gnd[k] - gga[k]; 
					//ff0[k] = fnd[k] - ffa[k]; gg0[k] = gnd[k] - gga[k];
				}
				if (mask[dde[k]] == 28 || mask[dde[k]] == 1)
				{
					ffh[k] = ffe[k]; ggh[k] = gge[k];
					ff0[k] += fne[k] - ffh[k]; gg0[k] += gne[k] - ggh[k];   
					//ff0[k] = fne[k] - ffh[k]; gg0[k] = gne[k] - ggh[k];
				}

				if (mask[ddf[k]] == 28 || mask[ddf[k]] == 1)
				{
					ffg[k] = fff[k]; ggg[k] = ggf[k];
					ff0[k] += fnf[k] - ffg[k]; gg0[k] += gnf[k] - ggg[k]; 
					//ff0[k] = fnf[k] - ffg[k]; gg0[k] = gnf[k] - ggg[k]; 
				}
				if (mask[ddg[k]] == 28 || mask[ddg[k]] == 1)
				{
					fff[k] = ffg[k]; ggf[k] = ggg[k];
					ff0[k] += fng[k] - fff[k]; gg0[k] += gng[k] - ggf[k];
					//ff0[k] = fng[k] - fff[k]; gg0[k] = gng[k] - ggf[k];
					
				}
				if (mask[ddh[k]] == 28 || mask[ddh[k]] == 1)
				{
					ffe[k] = ffh[k]; gge[k] = ggh[k];
					ff0[k] += fnh[k] - ffe[k]; gg0[k] += gnh[k] - gge[k]; 
					//ff0[k] = fnh[k] - ffe[k]; gg0[k] = gnh[k] - gge[k];  
				}
				if (mask[ddi[k]] == 28 || mask[ddi[k]] == 1)
				{
					ffl[k] = ffi[k]; ggl[k] = ggi[k];
					ff0[k] += fni[k] - ffl[k]; gg0[k] += gni[k] - ggl[k]; 
					//ff0[k] = fni[k] - ffl[k]; gg0[k] = gni[k] - ggl[k];
				}
				if (mask[ddj[k]] == 28 || mask[ddj[k]] == 1)
				{
					ffk[k] = ffj[k]; ggk[k] = ggj[k];
					ff0[k] += fnj[k] - ffk[k]; gg0[k] += gnj[k] - ggk[k]; 
					//ff0[k] = fnj[k] - ffk[k]; gg0[k] = gnj[k] - ggk[k];
				}
				if (mask[ddk[k]] == 28 || mask[ddk[k]] == 1)
				{
					ffj[k] = ffk[k]; ggj[k] = ggk[k];
					ff0[k] += fnk[k] - ffj[k]; gg0[k] += gnk[k] - ggj[k];  
					//ff0[k] = fnk[k] - ffj[k]; gg0[k] = gnk[k] - ggj[k];
				}	
				if (mask[ddl[k]] == 28 || mask[ddl[k]] == 1)
				{
					ffi[k] = ffl[k]; ggi[k] = ggl[k];
					ff0[k] += fnl[k] - ffi[k]; gg0[k] += gnl[k] - ggi[k];  
					//ff0[k] = fnl[k] - ffi[k]; gg0[k] = gnl[k] - ggi[k]; 
				}
				 
				break;

			case 28:
				break;	

			default:
			{
				computeCoordinate(k);
				cout << "In applyBoundaryConditions: A case is not defined at k = (" << xk << "," << yk << "," << zk << "), (" << (int) (mask[k]) << ")." << endl;
			}
				break;
		}
	}
}

//solid wall half way from a lattice row
/*		  ff2[k] = fn1[k]; gg2[k] = gn1[k];
 
 */	  
/*
if (mask[dd1[k]] == 28)
{       //solid wall is on th mask[k]=1 lattice nodes
	ff2[dd1[k]] = ff1[k]; gg2[dd1[k]] = gg1[k];
	//ff0[k] += fn1[k] - ff2[k]; gg0[k] += gn1[k] - gg2[k];   
	
	
}
if (mask[dd2[k]] == 28)
{
	ff1[dd2[k]] = ff2[k]; gg1[dd2[k]] = gg2[k];
	//ff0[k] += fn2[k] - ff1[k]; gg0[k] += gn2[k] - gg1[k];  
	//ff0[k] = fn2[k] - ff1[k]; gg0[k] = gn2[k] - gg1[k];
}
if (mask[dd3[k]] == 28)
{
	ff4[dd3[k]] = ff3[k];gg4[dd3[k]] = gg3[k];
	//ff0[k] += fn3[k] - ff4[k]; gg0[k] += gn3[k] - gg4[k];
	//ff0[k] = fn3[k] - ff4[k]; gg0[k] = gn3[k] - gg4[k];
}
if (mask[dd4[k]] == 28)
{
	ff3[dd4[k]] = ff4[k]; gg3[dd4[k]] = gg4[k];
	//ff0[k] += fn4[k] - ff3[k]; gg0[k] += gn4[k] - gg3[k];
	//ff0[k] = fn4[k] - ff3[k]; gg0[k] = gn4[k] - gg3[k];
}
if (mask[dd5[k]] == 28)
{
	ff6[dd5[k]] = ff5[k]; gg6[dd5[k]] = gg5[k];
	//ff0[k] += fn5[k] - ff6[k]; gg0[k] += gn5[k] - gg6[k];
	//ff0[k] = fn5[k] - ff6[k]; gg0[k] = gn5[k] - gg6[k];
}
if (mask[dd6[k]] == 28)
{
	ff5[dd6[k]] = ff6[k]; gg5[dd6[k]] = gg6[k];
	//ff0[k] += fn6[k] - ff5[k]; gg0[k] += gn6[k] - gg5[k];
	//ff0[k] = fn6[k] - ff5[k]; gg0[k] = gn6[k] - gg5[k];
}

if (mask[dda[k]] == 28)
{
	ffd[dda[k]] = ffa[k]; ggd[dda[k]] = gga[k];
	//ff0[k] += fna[k] - ffd[k]; gg0[k] += gna[k] - ggd[k]; 
	//ff0[k] = fna[k] - ffd[k]; gg0[k] = gna[k] - ggd[k];
}
if (mask[ddb[k]] == 28)
{
	ffc[ddb[k]] = ffb[k]; ggc[ddb[k]] = ggb[k];
	//ff0[k] += fnb[k] - ffc[k]; gg0[k] += gnb[k] - ggc[k];
	//ff0[k] = fnb[k] - ffc[k]; gg0[k] = gnb[k] - ggc[k];
}
if (mask[ddc[k]] == 28)
{
	ffb[ddc[k]] = ffc[k]; ggb[ddc[k]] = ggc[k];
	//ff0[k] += fnc[k] - ffb[k]; gg0[k] += gnc[k] - ggb[k]; 
	//ff0[k] = fnc[k] - ffb[k]; gg0[k] = gnc[k] - ggb[k];
}
if (mask[ddd[k]] == 28)
{	
	ffa[ddd[k]] = ffd[k]; gga[ddd[k]] = ggd[k];
	//ff0[k] += fnd[k] - ffa[k]; gg0[k] += gnd[k] - gga[k]; 
	//ff0[k] = fnd[k] - ffa[k]; gg0[k] = gnd[k] - gga[k];
}
if (mask[dde[k]] == 28)
{
	ffh[dde[k]] = ffe[k]; ggh[dde[k]] = gge[k];
	//ff0[k] += fne[k] - ffh[k]; gg0[k] += gne[k] - ggh[k];   
	//ff0[k] = fne[k] - ffh[k]; gg0[k] = gne[k] - ggh[k];
}

if (mask[ddf[k]] == 28)
{
	ffg[ddf[k]] = fff[k]; ggg[ddf[k]] = ggf[k];
	//ff0[k] += fnf[k] - ffg[k]; gg0[k] += gnf[k] - ggg[k]; 
	//ff0[k] = fnf[k] - ffg[k]; gg0[k] = gnf[k] - ggg[k]; 
}
if (mask[ddg[k]] == 28)
{
	fff[ddg[k]] = ffg[k]; ggf[ddg[k]] = ggg[k];
	//ff0[k] += fng[k] - fff[k]; gg0[k] += gng[k] - ggf[k];
	//ff0[k] = fng[k] - fff[k]; gg0[k] = gng[k] - ggf[k];
	
}
if (mask[ddh[k]] == 28)
{
	ffe[ddh[k]] = ffh[k]; gge[ddh[k]] = ggh[k];
	//ff0[k] += fnh[k] - ffe[k]; gg0[k] += gnh[k] - gge[k]; 
	//ff0[k] = fnh[k] - ffe[k]; gg0[k] = gnh[k] - gge[k];  
}
if (mask[ddi[k]] == 28)
{
	ffl[ddi[k]] = ffi[k]; ggl[ddi[k]] = ggi[k];
	//ff0[k] += fni[k] - ffl[k]; gg0[k] += gni[k] - ggl[k]; 
	//ff0[k] = fni[k] - ffl[k]; gg0[k] = gni[k] - ggl[k];
}
if (mask[ddj[k]] == 28)
{
	ffk[ddj[k]] = ffj[k]; ggk[ddj[k]] = ggj[k];
	//ff0[k] += fnj[k] - ffk[k]; gg0[k] += gnj[k] - ggk[k]; 
	//ff0[k] = fnj[k] - ffk[k]; gg0[k] = gnj[k] - ggk[k];
}
if (mask[ddk[k]] == 28)
{
	ffj[ddk[k]] = ffk[k]; ggj[ddk[k]] = ggk[k];
	//ff0[k] += fnk[k] - ffj[k]; gg0[k] += gnk[k] - ggj[k];  
	//ff0[k] = fnk[k] - ffj[k]; gg0[k] = gnk[k] - ggj[k];
}	
if (mask[ddl[k]] == 28)
{
	ffi[ddl[k]] = ffl[k]; ggi[ddl[k]] = ggl[k];
	//ff0[k] += fnl[k] - ffi[k]; gg0[k] += gnl[k] - ggi[k];  
	//ff0[k] = fnl[k] - ffi[k]; gg0[k] = gnl[k] - ggi[k]; 
}
*/

