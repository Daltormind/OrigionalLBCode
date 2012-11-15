//densitiesAtSolidBoundaries.cpp: initialise the density of fluid g at the boundaries
//Marco, 30 June 2011
//Matthew This sets up the appropriate boundry values of the phase to satisfy the wetting boundry condition

#include "wet.h"

void wet::densitiesAtSolidBoundaries(void)
{

	//if(t%infoStep==0)
	//	cout << "Process " << rank << ": setting densities at solid boundaries..." << endl;

	double sqrt2 = sqrt(2.0L);
	double sqrt3 = sqrt(3.0L);	
	double xi = sqrt(8.0L*kappa_p/A);
	
	double pEq = 0.0; // all these pOld,pEq,pNew are in the 28 layer
	double pOld = 0.0;
	double pNew = 0.0;
	
	for( k = k1 ; k < k2 ; k++)
	{
		if (mask[k] == 1) 
		{	
			
			if (mask[dd1[k]] == 28){
				
				//				_____
				//right wall-->|	 |
				//			___|     |___
	
				dp_dx = ((p[dd1[k]] - p[dd2[k]])/6 + (p[dda[k]]+p[ddc[k]]+p[ddi[k]]+p[ddk[k]]-p[ddb[k]]-p[ddd[k]]-p[ddj[k]]-p[ddl[k]])/12);
				//dp_dy = ((p[dd3[k]] - p[dd4[k]])/6 + (p[dda[k]]+p[ddb[k]]+p[dde[k]]+p[ddg[k]]-p[ddc[k]]-p[ddd[k]]-p[ddf[k]]-p[ddh[k]])/12);
				//dp_dz = ((p[dd5[k]] - p[dd6[k]])/6 + (p[dde[k]]+p[ddf[k]]+p[ddi[k]]+p[ddj[k]]-p[ddg[k]]-p[ddh[k]]-p[ddk[k]]-p[ddl[k]])/12);
				
				pEq = p[dd2[k]] - 2*phi11;   
				pOld = p[dd2[k]] + 2*dp_dx;  
				pNew = pOld - (pOld - pEq)/tausurf;
				p[dd1[k]] = pNew;
				
			}
			
			if (mask[dd2[k]] == 28) 
			{
				//				_____
				//left wall    |	 |<--
				//			___|     |___
				
				dp_dx = ((p[dd1[k]] - p[dd2[k]])/6 + (p[dda[k]]+p[ddc[k]]+p[ddi[k]]+p[ddk[k]]-p[ddb[k]]-p[ddd[k]]-p[ddj[k]]-p[ddl[k]])/12);
				//dp_dy = ((p[dd3[k]] - p[dd4[k]])/6 + (p[dda[k]]+p[ddb[k]]+p[dde[k]]+p[ddg[k]]-p[ddc[k]]-p[ddd[k]]-p[ddf[k]]-p[ddh[k]])/12);
				//dp_dz = ((p[dd5[k]] - p[dd6[k]])/6 + (p[dde[k]]+p[ddf[k]]+p[ddi[k]]+p[ddj[k]]-p[ddg[k]]-p[ddh[k]]-p[ddk[k]]-p[ddl[k]])/12);
				
				pEq = p[dd1[k]] - 2*phi11;   
				pOld = p[dd1[k]] - 2*dp_dx;  				
				pNew = pOld - (pOld - pEq)/tausurf;
				p[dd2[k]] = pNew;
				
			}
 			if (mask[dd6[k]] == 28)			
			{										//				__v__
				//p[dd6[k]] = p[dd5[k]] - 2*phi11;	//bottom wall  |     |
													//			___|     |___
			
				//dp_dx = ((p[dd1[k]] - p[dd2[k]])/6 + (p[dda[k]]+p[ddc[k]]+p[ddi[k]]+p[ddk[k]]-p[ddb[k]]-p[ddd[k]]-p[ddj[k]]-p[ddl[k]])/12);
				//dp_dy = ((p[dd3[k]] - p[dd4[k]])/6 + (p[dda[k]]+p[ddb[k]]+p[dde[k]]+p[ddg[k]]-p[ddc[k]]-p[ddd[k]]-p[ddf[k]]-p[ddh[k]])/12);
				dp_dz = ((p[dd5[k]] - p[dd6[k]])/6 + (p[dde[k]]+p[ddf[k]]+p[ddi[k]]+p[ddj[k]]-p[ddg[k]]-p[ddh[k]]-p[ddk[k]]-p[ddl[k]])/12);
				
				pEq = p[dd5[k]] - 2*phi11;   
				pOld = p[dd5[k]] - 2*dp_dz;  
				pNew = pOld - (pOld - pEq)/tausurf;
				p[dd6[k]] = pNew;
				
			}
			if(mask[dd5[k]] == 28)
			{										//			---|     |---
				//p[dd5[k]] = p[dd6[k]] - 2*phi11;	//above wall   |_____|
													//				  A	
				
				//dp_dx = ((p[dd1[k]] - p[dd2[k]])/6 + (p[dda[k]]+p[ddc[k]]+p[ddi[k]]+p[ddk[k]]-p[ddb[k]]-p[ddd[k]]-p[ddj[k]]-p[ddl[k]])/12);
				//dp_dy = ((p[dd3[k]] - p[dd4[k]])/6 + (p[dda[k]]+p[ddb[k]]+p[dde[k]]+p[ddg[k]]-p[ddc[k]]-p[ddd[k]]-p[ddf[k]]-p[ddh[k]])/12);
				dp_dz = ((p[dd5[k]] - p[dd6[k]])/6 + (p[dde[k]]+p[ddf[k]]+p[ddi[k]]+p[ddj[k]]-p[ddg[k]]-p[ddh[k]]-p[ddk[k]]-p[ddl[k]])/12);
				
				pEq = p[dd6[k]] - 2*phi11;   // all these pOld,pEq,pNew are in the 28 layer
				pOld = p[dd6[k]] + 2*dp_dz;  
				pNew = pOld - (pOld - pEq)/tausurf;
				p[dd5[k]] = pNew;
				
			}
		
			
			if (mask[dd3[k]] == 28) 
			{										//				_____
				//p[dd3[k]] = p[dd4[k]] - 2*phi11;	//front wall   |     |	
													//			___|  x  |___
				
				//dp_dx = ((p[dd1[k]] - p[dd2[k]])/6 + (p[dda[k]]+p[ddc[k]]+p[ddi[k]]+p[ddk[k]]-p[ddb[k]]-p[ddd[k]]-p[ddj[k]]-p[ddl[k]])/12);
				dp_dy = ((p[dd3[k]] - p[dd4[k]])/6 + (p[dda[k]]+p[ddb[k]]+p[dde[k]]+p[ddg[k]]-p[ddc[k]]-p[ddd[k]]-p[ddf[k]]-p[ddh[k]])/12);
				//dp_dz = ((p[dd5[k]] - p[dd6[k]])/6 + (p[dde[k]]+p[ddf[k]]+p[ddi[k]]+p[ddj[k]]-p[ddg[k]]-p[ddh[k]]-p[ddk[k]]-p[ddl[k]])/12);
				
				pEq = p[dd4[k]] - 2*phi11;   // all these pOld,pEq,pNew are in the 28 layer
				pOld = p[dd4[k]] + 2*dp_dy;  
				pNew = pOld - (pOld - pEq)/tausurf;
				p[dd3[k]] = pNew;
				
			}
			if (mask[dd4[k]] == 28) 
			{										//				_____
				//p[dd4[k]] = p[dd3[k]] - 2*phi11;	//back wall	   |     |
													//			___|  o  |___
				
				//dp_dx = ((p[dd1[k]] - p[dd2[k]])/6 + (p[dda[k]]+p[ddc[k]]+p[ddi[k]]+p[ddk[k]]-p[ddb[k]]-p[ddd[k]]-p[ddj[k]]-p[ddl[k]])/12);
				dp_dy = ((p[dd3[k]] - p[dd4[k]])/6 + (p[dda[k]]+p[ddb[k]]+p[dde[k]]+p[ddg[k]]-p[ddc[k]]-p[ddd[k]]-p[ddf[k]]-p[ddh[k]])/12);
				//dp_dz = ((p[dd5[k]] - p[dd6[k]])/6 + (p[dde[k]]+p[ddf[k]]+p[ddi[k]]+p[ddj[k]]-p[ddg[k]]-p[ddh[k]]-p[ddk[k]]-p[ddl[k]])/12);
				
				pEq = p[dd3[k]] - 2*phi11;   // all these pOld,pEq,pNew are in the 28 layer
				pOld = p[dd3[k]] - 2*dp_dy;  
				pNew = pOld - (pOld - pEq)/tausurf;
				p[dd4[k]] = pNew;
				
			}
			
			//old:
			/*
			if (mask[dd1[k]] == 28)
			{										//				_____
				p[dd1[k]] = p[dd2[k]] - 2*phi11;	//right wall-->|	 |
			}										//			___|     |___
			if (mask[dd2[k]] == 28) 
			{										//				_____
				p[dd2[k]] = p[dd1[k]] - 2*phi11;	//left wall	   |     |<--
			}										//			___|     |___
 			if (mask[dd6[k]] == 28)			
			{										//				__v__
				p[dd6[k]] = p[dd5[k]] - 2*phi11;	//bottom wall  |     |
			}										//			___|     |___
			if(mask[dd5[k]] == 28)
			{										//			---|     |---
				p[dd5[k]] = p[dd6[k]] - 2*phi11;	//above wall   |_____|
			}										//				  A	
			if (mask[dd3[k]] == 28) 
			{										//				_____
				p[dd3[k]] = p[dd4[k]] - 2*phi11;	//front wall   |     |	
			}										//			___|  x  |___
			if (mask[dd4[k]] == 28) 
			{										//				_____
				p[dd4[k]] = p[dd3[k]] - 2*phi11;	//back wall	   |     |
			}										//			___|  o  |___
			*/                                        		     
			
		}//endif(mask[k]==1)
	
	}
	
	//repeat so that edges are written AFTER walls
	double factor = 0.5;
	double factor2 = 0;
	for( k = k1 ; k < k2 ; k++)
	{
	
		if (mask[k]==1) {
			
			if(mask[dd1[k]]==1 && mask[dd6[k]]==1 && mask[ddk[k]]==28)
			{	
				//				__v__
				//edge A	-->|     |
				//	 		___|     |__
				
				//p[ddk[k]] = factor*(p[ddi[k]] - 2*phi11 + p[ddl[k]] - 2*phi11) + factor2*(p[ddj[k]] - 2*sqrt2*phi11);
				
				dp_dx = ((p[dd1[k]] - p[dd2[k]])/6 + (p[dda[k]]+p[ddc[k]]+p[ddi[k]]+p[ddk[k]]-p[ddb[k]]-p[ddd[k]]-p[ddj[k]]-p[ddl[k]])/12);
				//dp_dy = ((p[dd3[k]] - p[dd4[k]])/6 + (p[dda[k]]+p[ddb[k]]+p[dde[k]]+p[ddg[k]]-p[ddc[k]]-p[ddd[k]]-p[ddf[k]]-p[ddh[k]])/12);
				dp_dz = ((p[dd5[k]] - p[dd6[k]])/6 + (p[dde[k]]+p[ddf[k]]+p[ddi[k]]+p[ddj[k]]-p[ddg[k]]-p[ddh[k]]-p[ddk[k]]-p[ddl[k]])/12);
				
				pEq =  factor*(p[ddi[k]] - 2*phi11 + p[ddl[k]] - 2*phi11) + factor2*(p[ddj[k]] - 2*sqrt2*phi11);
				pOld = p[ddj[k]] + dp_dx - dp_dz;  
				pNew = pOld - (pOld - pEq)/tausurf;
				p[ddk[k]] = pNew;
				
			}	
			if(mask[dd2[k]]==1 && mask[dd6[k]]==1 && mask[ddl[k]]==28)
			{	
				//				__v__	
				//edge B	   |     |<--
				//			___|     |___	
				
				//p[ddl[k]] = factor*(p[ddj[k]] - 2*phi11 + p[ddk[k]] - 2*phi11) + factor2*(p[ddi[k]] - 2*sqrt2*phi11);
				
				dp_dx = ((p[dd1[k]] - p[dd2[k]])/6 + (p[dda[k]]+p[ddc[k]]+p[ddi[k]]+p[ddk[k]]-p[ddb[k]]-p[ddd[k]]-p[ddj[k]]-p[ddl[k]])/12);
				//dp_dy = ((p[dd3[k]] - p[dd4[k]])/6 + (p[dda[k]]+p[ddb[k]]+p[dde[k]]+p[ddg[k]]-p[ddc[k]]-p[ddd[k]]-p[ddf[k]]-p[ddh[k]])/12);
				dp_dz = ((p[dd5[k]] - p[dd6[k]])/6 + (p[dde[k]]+p[ddf[k]]+p[ddi[k]]+p[ddj[k]]-p[ddg[k]]-p[ddh[k]]-p[ddk[k]]-p[ddl[k]])/12);
				
				pEq =  factor*(p[ddj[k]] - 2*phi11 + p[ddk[k]] - 2*phi11) + factor2*(p[ddi[k]] - 2*sqrt2*phi11);
				pOld = p[ddi[k]] - dp_dx - dp_dz;  
				pNew = pOld - (pOld - pEq)/tausurf;
				p[ddl[k]] = pNew;
				
				
			}
			if(mask[dd3[k]]==1 && mask[dd6[k]]==1 && mask[ddg[k]]==28)
			{		
				//				__v__
				//edge C	   |     |
				//			___|  x  |___
				
				//p[ddg[k]] = factor*(p[dde[k]] - 2*phi11 + p[ddh[k]] - 2*phi11) + factor2*(p[ddf[k]] - 2*sqrt2*phi11); 
				
				//dp_dx = ((p[dd1[k]] - p[dd2[k]])/6 + (p[dda[k]]+p[ddc[k]]+p[ddi[k]]+p[ddk[k]]-p[ddb[k]]-p[ddd[k]]-p[ddj[k]]-p[ddl[k]])/12);
				dp_dy = ((p[dd3[k]] - p[dd4[k]])/6 + (p[dda[k]]+p[ddb[k]]+p[dde[k]]+p[ddg[k]]-p[ddc[k]]-p[ddd[k]]-p[ddf[k]]-p[ddh[k]])/12);
				dp_dz = ((p[dd5[k]] - p[dd6[k]])/6 + (p[dde[k]]+p[ddf[k]]+p[ddi[k]]+p[ddj[k]]-p[ddg[k]]-p[ddh[k]]-p[ddk[k]]-p[ddl[k]])/12);
				
				pEq =  factor*(p[dde[k]] - 2*phi11 + p[ddh[k]] - 2*phi11) + factor2*(p[ddf[k]] - 2*sqrt2*phi11); 
				pOld = p[ddf[k]] + dp_dy - dp_dz;  
				pNew = pOld - (pOld - pEq)/tausurf;
				p[ddg[k]] = pNew;
				
			}
			if(mask[dd4[k]]==1 && mask[dd6[k]]==1 && mask[ddh[k]]==28)
			{	
				//				__v__
				//edge D	   |     |
				//			___|  o  |___
				
				//p[ddh[k]] = factor*(p[ddf[k]] - 2*phi11 + p[ddg[k]] - 2*phi11) + factor2*(p[dde[k]] - 2*sqrt2*phi11); 
				
				//dp_dx = ((p[dd1[k]] - p[dd2[k]])/6 + (p[dda[k]]+p[ddc[k]]+p[ddi[k]]+p[ddk[k]]-p[ddb[k]]-p[ddd[k]]-p[ddj[k]]-p[ddl[k]])/12);
				dp_dy = ((p[dd3[k]] - p[dd4[k]])/6 + (p[dda[k]]+p[ddb[k]]+p[dde[k]]+p[ddg[k]]-p[ddc[k]]-p[ddd[k]]-p[ddf[k]]-p[ddh[k]])/12);
				dp_dz = ((p[dd5[k]] - p[dd6[k]])/6 + (p[dde[k]]+p[ddf[k]]+p[ddi[k]]+p[ddj[k]]-p[ddg[k]]-p[ddh[k]]-p[ddk[k]]-p[ddl[k]])/12);
				
				pEq =  factor*(p[ddf[k]] - 2*phi11 + p[ddg[k]] - 2*phi11) + factor2*(p[dde[k]] - 2*sqrt2*phi11); 
				pOld = p[dde[k]] - dp_dy - dp_dz;  
				pNew = pOld - (pOld - pEq)/tausurf;
				p[ddh[k]] = pNew;
				
			}
			
			//  around post edges
			if(mask[dd1[k]]==1 && mask[dd3[k]]==1 && mask[dda[k]]==28)
			{											
				//				_____
				//edge E	-->|     |
				//			___|  x  |___
				
				//p[dda[k]] = factor*(p[ddb[k]] - 2*phi11 + p[ddc[k]] - 2*phi11) + factor2*(p[ddd[k]] - 2*sqrt2*phi11);
				
				dp_dx = ((p[dd1[k]] - p[dd2[k]])/6 + (p[dda[k]]+p[ddc[k]]+p[ddi[k]]+p[ddk[k]]-p[ddb[k]]-p[ddd[k]]-p[ddj[k]]-p[ddl[k]])/12);
				dp_dy = ((p[dd3[k]] - p[dd4[k]])/6 + (p[dda[k]]+p[ddb[k]]+p[dde[k]]+p[ddg[k]]-p[ddc[k]]-p[ddd[k]]-p[ddf[k]]-p[ddh[k]])/12);
				//dp_dz = ((p[dd5[k]] - p[dd6[k]])/6 + (p[dde[k]]+p[ddf[k]]+p[ddi[k]]+p[ddj[k]]-p[ddg[k]]-p[ddh[k]]-p[ddk[k]]-p[ddl[k]])/12);
				
				pEq =  factor*(p[ddb[k]] - 2*phi11 + p[ddc[k]] - 2*phi11) + factor2*(p[ddd[k]] - 2*sqrt2*phi11);
				pOld = p[ddd[k]] + dp_dx + dp_dy;  
				pNew = pOld - (pOld - pEq)/tausurf;
				p[dda[k]] = pNew;
				
				
			}
			if(mask[dd2[k]]==1 && mask[dd3[k]]==1 && mask[ddb[k]]==28)
			{											
				//				_____
				//edge F	   |     |<--
				//			___|  x  |___
				
				//p[ddb[k]] = factor*(p[ddd[k]] - 2*phi11 + p[dda[k]] - 2*phi11) + factor2*(p[ddc[k]] - 2*sqrt2*phi11);
				
				dp_dx = ((p[dd1[k]] - p[dd2[k]])/6 + (p[dda[k]]+p[ddc[k]]+p[ddi[k]]+p[ddk[k]]-p[ddb[k]]-p[ddd[k]]-p[ddj[k]]-p[ddl[k]])/12);
				dp_dy = ((p[dd3[k]] - p[dd4[k]])/6 + (p[dda[k]]+p[ddb[k]]+p[dde[k]]+p[ddg[k]]-p[ddc[k]]-p[ddd[k]]-p[ddf[k]]-p[ddh[k]])/12);
				//dp_dz = ((p[dd5[k]] - p[dd6[k]])/6 + (p[dde[k]]+p[ddf[k]]+p[ddi[k]]+p[ddj[k]]-p[ddg[k]]-p[ddh[k]]-p[ddk[k]]-p[ddl[k]])/12);
				
				pEq =  factor*(p[ddd[k]] - 2*phi11 + p[dda[k]] - 2*phi11) + factor2*(p[ddc[k]] - 2*sqrt2*phi11);
				pOld = p[ddc[k]] - dp_dx + dp_dy;  
				pNew = pOld - (pOld - pEq)/tausurf;
				p[ddb[k]] = pNew;
				
			}	
			if(mask[dd4[k]]==1 && mask[dd2[k]]==1 && mask[ddd[k]]==28)
			{	
				//				_____
				//edge G	   |     |<--	
				//			___|  o  |___
				
				//p[ddd[k]] = factor*(p[ddb[k]] - 2*phi11 + p[ddc[k]] - 2*phi11) + factor2*(p[dda[k]] - 2*sqrt2*phi11);
				
				dp_dx = ((p[dd1[k]] - p[dd2[k]])/6 + (p[dda[k]]+p[ddc[k]]+p[ddi[k]]+p[ddk[k]]-p[ddb[k]]-p[ddd[k]]-p[ddj[k]]-p[ddl[k]])/12);
				dp_dy = ((p[dd3[k]] - p[dd4[k]])/6 + (p[dda[k]]+p[ddb[k]]+p[dde[k]]+p[ddg[k]]-p[ddc[k]]-p[ddd[k]]-p[ddf[k]]-p[ddh[k]])/12);
				//dp_dz = ((p[dd5[k]] - p[dd6[k]])/6 + (p[dde[k]]+p[ddf[k]]+p[ddi[k]]+p[ddj[k]]-p[ddg[k]]-p[ddh[k]]-p[ddk[k]]-p[ddl[k]])/12);
				
				pEq =  factor*(p[ddb[k]] - 2*phi11 + p[ddc[k]] - 2*phi11) + factor2*(p[dda[k]] - 2*sqrt2*phi11);
				pOld = p[dda[k]] - dp_dx - dp_dy;  
				pNew = pOld - (pOld - pEq)/tausurf;
				p[ddd[k]] = pNew;
				
			}
			if(mask[dd1[k]]==1 && mask[dd4[k]]==1 && mask[ddc[k]]==28)
			{	
				//				_____
				//edge H	-->|     |
				//			___|  o  |___
				
				//p[ddc[k]] = factor*(p[ddd[k]] - 2*phi11 + p[dda[k]] - 2*phi11) + factor2*(p[ddb[k]] - 2*sqrt2*phi11);
				
				dp_dx = ((p[dd1[k]] - p[dd2[k]])/6 + (p[dda[k]]+p[ddc[k]]+p[ddi[k]]+p[ddk[k]]-p[ddb[k]]-p[ddd[k]]-p[ddj[k]]-p[ddl[k]])/12);
				dp_dy = ((p[dd3[k]] - p[dd4[k]])/6 + (p[dda[k]]+p[ddb[k]]+p[dde[k]]+p[ddg[k]]-p[ddc[k]]-p[ddd[k]]-p[ddf[k]]-p[ddh[k]])/12);
				//dp_dz = ((p[dd5[k]] - p[dd6[k]])/6 + (p[dde[k]]+p[ddf[k]]+p[ddi[k]]+p[ddj[k]]-p[ddg[k]]-p[ddh[k]]-p[ddk[k]]-p[ddl[k]])/12);
				
				pEq =  factor*(p[ddd[k]] - 2*phi11 + p[dda[k]] - 2*phi11) + factor2*(p[ddb[k]] - 2*sqrt2*phi11);
				pOld = p[ddb[k]] + dp_dx - dp_dy;  
				pNew = pOld - (pOld - pEq)/tausurf;
				p[ddc[k]] = pNew;
				
			}
			
			// below posts edges
			if(mask[dd1[k]]==1 && mask[dd5[k]]==1 && mask[ddi[k]]==28)
			{	
				//			---|     |---
				//edge I	-->|_____|
				//		        A	
				
				//p[ddi[k]] = factor*(p[ddj[k]] - 2*phi11 + p[ddk[k]] - 2*phi11) + factor2*(p[ddl[k]] - 2*sqrt2*phi11); 
				
				dp_dx = ((p[dd1[k]] - p[dd2[k]])/6 + (p[dda[k]]+p[ddc[k]]+p[ddi[k]]+p[ddk[k]]-p[ddb[k]]-p[ddd[k]]-p[ddj[k]]-p[ddl[k]])/12);
				//dp_dy = ((p[dd3[k]] - p[dd4[k]])/6 + (p[dda[k]]+p[ddb[k]]+p[dde[k]]+p[ddg[k]]-p[ddc[k]]-p[ddd[k]]-p[ddf[k]]-p[ddh[k]])/12);
				dp_dz = ((p[dd5[k]] - p[dd6[k]])/6 + (p[dde[k]]+p[ddf[k]]+p[ddi[k]]+p[ddj[k]]-p[ddg[k]]-p[ddh[k]]-p[ddk[k]]-p[ddl[k]])/12);
				
				pEq =  factor*(p[ddj[k]] - 2*phi11 + p[ddk[k]] - 2*phi11) + factor2*(p[ddl[k]] - 2*sqrt2*phi11); 
				pOld = p[ddl[k]] + dp_dx + dp_dz;  
				pNew = pOld - (pOld - pEq)/tausurf;
				p[ddi[k]] = pNew;
				
			}
			if(mask[dd2[k]]==1 && mask[dd5[k]]==1 && mask[ddj[k]]==28)
			{	
				//			---|     |---
				//edge M	   |_____|<--
				//		            A
				
				//p[ddj[k]] = factor*(p[ddi[k]] - 2*phi11 + p[ddl[k]] - 2*phi11) + factor2*(p[ddk[k]] - 2*sqrt2*phi11); 
				
				dp_dx = ((p[dd1[k]] - p[dd2[k]])/6 + (p[dda[k]]+p[ddc[k]]+p[ddi[k]]+p[ddk[k]]-p[ddb[k]]-p[ddd[k]]-p[ddj[k]]-p[ddl[k]])/12);
				//dp_dy = ((p[dd3[k]] - p[dd4[k]])/6 + (p[dda[k]]+p[ddb[k]]+p[dde[k]]+p[ddg[k]]-p[ddc[k]]-p[ddd[k]]-p[ddf[k]]-p[ddh[k]])/12);
				dp_dz = ((p[dd5[k]] - p[dd6[k]])/6 + (p[dde[k]]+p[ddf[k]]+p[ddi[k]]+p[ddj[k]]-p[ddg[k]]-p[ddh[k]]-p[ddk[k]]-p[ddl[k]])/12);
				
				pEq =  factor*(p[ddi[k]] - 2*phi11 + p[ddl[k]] - 2*phi11) + factor2*(p[ddk[k]] - 2*sqrt2*phi11);
				pOld = p[ddk[k]] - dp_dx + dp_dz;  
				pNew = pOld - (pOld - pEq)/tausurf;
				p[ddj[k]] = pNew;
			
			}
			if(mask[dd4[k]]==1 && mask[dd5[k]]==1 && mask[ddf[k]]==28)
			{	
				//			---|  o  |---
				//edge N	   |_____|
				//	   	          A
				
				//p[ddf[k]] = factor*(p[dde[k]] - 2*phi11 + p[ddh[k]] - 2*phi11) + factor2*(p[ddg[k]] - 2*sqrt2*phi11); 
				
				//dp_dx = ((p[dd1[k]] - p[dd2[k]])/6 + (p[dda[k]]+p[ddc[k]]+p[ddi[k]]+p[ddk[k]]-p[ddb[k]]-p[ddd[k]]-p[ddj[k]]-p[ddl[k]])/12);
				dp_dy = ((p[dd3[k]] - p[dd4[k]])/6 + (p[dda[k]]+p[ddb[k]]+p[dde[k]]+p[ddg[k]]-p[ddc[k]]-p[ddd[k]]-p[ddf[k]]-p[ddh[k]])/12);
				dp_dz = ((p[dd5[k]] - p[dd6[k]])/6 + (p[dde[k]]+p[ddf[k]]+p[ddi[k]]+p[ddj[k]]-p[ddg[k]]-p[ddh[k]]-p[ddk[k]]-p[ddl[k]])/12);
				
				pEq =  factor*(p[dde[k]] - 2*phi11 + p[ddh[k]] - 2*phi11) + factor2*(p[ddg[k]] - 2*sqrt2*phi11);
				pOld = p[ddg[k]] - dp_dy + dp_dz;  
				pNew = pOld - (pOld - pEq)/tausurf;
				p[ddf[k]] = pNew;
				
			}
			if(mask[dd3[k]]==1 && mask[dd5[k]]==1 && mask[dde[k]]==28)
			{	
				//			---|  x  |---
				//edge  	   |_____|
				//		          A
				
				//p[dde[k]] = factor*(p[ddf[k]] - 2*phi11 + p[ddg[k]] - 2*phi11) + factor2*(p[ddh[k]] - 2*sqrt2*phi11);
				
				//dp_dx = ((p[dd1[k]] - p[dd2[k]])/6 + (p[dda[k]]+p[ddc[k]]+p[ddi[k]]+p[ddk[k]]-p[ddb[k]]-p[ddd[k]]-p[ddj[k]]-p[ddl[k]])/12);
				dp_dy = ((p[dd3[k]] - p[dd4[k]])/6 + (p[dda[k]]+p[ddb[k]]+p[dde[k]]+p[ddg[k]]-p[ddc[k]]-p[ddd[k]]-p[ddf[k]]-p[ddh[k]])/12);
				dp_dz = ((p[dd5[k]] - p[dd6[k]])/6 + (p[dde[k]]+p[ddf[k]]+p[ddi[k]]+p[ddj[k]]-p[ddg[k]]-p[ddh[k]]-p[ddk[k]]-p[ddl[k]])/12);
				
				pEq =  factor*(p[ddf[k]] - 2*phi11 + p[ddg[k]] - 2*phi11) + factor2*(p[ddh[k]] - 2*sqrt2*phi11);
				pOld = p[ddh[k]] + dp_dy + dp_dz;  
				pNew = pOld - (pOld - pEq)/tausurf;
				p[dde[k]] = pNew;
				
				
			}
		 /*
		 // **** OLD
			//  top post edges
			if(mask[dd1[k]]==1 && mask[dd6[k]]==1 && mask[ddk[k]]==28)
			{	
				//				__v__
				//edge A	-->|     |
				//	 		___|     |__
				
				p[ddk[k]] = factor*(p[ddi[k]] - 2*phi11 + p[ddl[k]] - 2*phi11) + factor2*(p[ddj[k]] - 2*sqrt2*phi11);
			}	
			if(mask[dd2[k]]==1 && mask[dd6[k]]==1 && mask[ddl[k]]==28)
			{	
				//				__v__	
				//edge B	   |     |<--
				//			___|     |___	
				
				p[ddl[k]] = factor*(p[ddj[k]] - 2*phi11 + p[ddk[k]] - 2*phi11) + factor2*(p[ddi[k]] - 2*sqrt2*phi11);
			}
			if(mask[dd3[k]]==1 && mask[dd6[k]]==1 && mask[ddg[k]]==28)
			{		
				//				__v__
				//edge C	   |     |
				//			___|  x  |___
				
				p[ddg[k]] = factor*(p[dde[k]] - 2*phi11 + p[ddh[k]] - 2*phi11) + factor2*(p[ddf[k]] - 2*sqrt2*phi11); 
			}
			if(mask[dd4[k]]==1 && mask[dd6[k]]==1 && mask[ddh[k]]==28)
			{	
				//				__v__
				//edge D	   |     |
				//			___|  o  |___
				
				p[ddh[k]] = factor*(p[ddf[k]] - 2*phi11 + p[ddg[k]] - 2*phi11) + factor2*(p[dde[k]] - 2*sqrt2*phi11); 
			}
			
			//  around post edges
			if(mask[dd1[k]]==1 && mask[dd3[k]]==1 && mask[dda[k]]==28)
			{											
				//				_____
				//edge E	-->|     |
				//			___|  x  |___
				
				p[dda[k]] = factor*(p[ddb[k]] - 2*phi11 + p[ddc[k]] - 2*phi11) + factor2*(p[ddd[k]] - 2*sqrt2*phi11);
			}
			if(mask[dd2[k]]==1 && mask[dd3[k]]==1 && mask[ddb[k]]==28)
			{											
				//				_____
				//edge F	   |     |<--
				//			___|  x  |___
				
				p[ddb[k]] = factor*(p[ddd[k]] - 2*phi11 + p[dda[k]] - 2*phi11) + factor2*(p[ddc[k]] - 2*sqrt2*phi11);
			}	
			if(mask[dd4[k]]==1 && mask[dd2[k]]==1 && mask[ddd[k]]==28)
			{	
				//				_____
				//edge G	   |     |<--	
				//			___|  o  |___
				
				p[ddd[k]] = factor*(p[ddb[k]] - 2*phi11 + p[ddc[k]] - 2*phi11) + factor2*(p[dda[k]] - 2*sqrt2*phi11);
			}
			if(mask[dd1[k]]==1 && mask[dd4[k]]==1 && mask[ddc[k]]==28)
			{	
				//				_____
				//edge H	-->|     |
				//			___|  o  |___
				
				p[ddc[k]] = factor*(p[ddd[k]] - 2*phi11 + p[dda[k]] - 2*phi11) + factor2*(p[ddb[k]] - 2*sqrt2*phi11);
			}
			
			// below posts edges
			if(mask[dd1[k]]==1 && mask[dd5[k]]==1 && mask[ddi[k]]==28)
			{	
				//			---|     |---
				//edge I	-->|_____|
				//		        A	
				p[ddi[k]] = factor*(p[ddj[k]] - 2*phi11 + p[ddk[k]] - 2*phi11) + factor2*(p[ddl[k]] - 2*sqrt2*phi11); 
			}
			if(mask[dd2[k]]==1 && mask[dd5[k]]==1 && mask[ddj[k]]==28)
			{	
				//			---|     |---
				//edge M	   |_____|<--
				//		            A
				
				p[ddj[k]] = factor*(p[ddi[k]] - 2*phi11 + p[ddl[k]] - 2*phi11) + factor2*(p[ddk[k]] - 2*sqrt2*phi11); 
				
			}
			if(mask[dd4[k]]==1 && mask[dd5[k]]==1 && mask[ddf[k]]==28)
			{	
				//			---|  o  |---
				//edge N	   |_____|
				//	   	          A
				
				p[ddf[k]] = factor*(p[dde[k]] - 2*phi11 + p[ddh[k]] - 2*phi11) + factor2*(p[ddg[k]] - 2*sqrt2*phi11); 
			}
			if(mask[dd3[k]]==1 && mask[dd5[k]]==1 && mask[dde[k]]==28)
			{	
				//			---|  x  |---
				//edge  	   |_____|
				//		          A
				
				p[dde[k]] = factor*(p[ddf[k]] - 2*phi11 + p[ddg[k]] - 2*phi11) + factor2*(p[ddh[k]] - 2*sqrt2*phi11);
			}
		  */
		}	
	}	
	//if(t%infoStep==0)
	//		cout << "Process " << rank << ": densities at solid boundaries set." << endl;
}



/*
 // 4 top vertices	
 if(mask[dd1[k]]==1 && mask[dd4[k]]==1 && mask[dd6[k]]==1 && mask[dd1[dd4[dd6[k]]]]==28)	
 {								//	    __v__
 p[dd1[dd4[dd6[k]]]] = p[dd2[dd3[dd5[k]]]] - 2*sqrt3*phi11;//	-->|     |
 }								//	___|  o  |___
 if(mask[dd1[k]]==1 && mask[dd3[k]]==1 && mask[dd6[k]]==1 && mask[dd1[dd3[dd6[k]]]]==28)
 {								//	    __v__
 p[dd1[dd3[dd6[k]]]] = p[dd2[dd4[dd5[k]]]] - 2*sqrt3*phi11;//	-->|     |
 }								//	___|  x  |___
 if(mask[dd2[k]]==1 && mask[dd3[k]]==1 && mask[dd6[k]]==1 && mask[dd2[dd3[dd6[k]]]]==28) 
 {								//	    __v__
 p[dd2[dd3[dd6[k]]]] = p[dd1[dd4[dd5[k]]]] - 2*sqrt3*phi11;//	   |     |<--
 }								//	___|  x  |___
 if(mask[dd2[k]]==1 && mask[dd4[k]]==1 && mask[dd6[k]]==1 && mask[dd2[dd4[dd6[k]]]]==28)	
 {								//	    __v__
 p[dd2[dd4[dd6[k]]]] = p[dd1[dd3[dd5[k]]]] - 2*sqrt3*phi11;//	   |     |<--
 }								//	___|  o  |___
 
 
 
 
 // 4 bottom vertices
 if(mask[dd1[k]]==1 && mask[dd4[k]]==1 && mask[dd5[k]]==1 && mask[dd1[dd4[dd5[k]]]]==28)	
 {							//		---|  o  |---
 p[dd1[dd4[dd5[k]]]] = p[dd2[dd3[dd6[k]]]] - 2*sqrt3*phi11;//	-->|_____|
 }							//		      A
 if(mask[dd1[k]]==1 && mask[dd3[k]]==1 && mask[dd5[k]]==1 && mask[dd1[dd3[dd5[k]]]]==28)	
 {							//		---|  x  |---
 p[dd1[dd3[dd5[k]]]] = p[dd2[dd4[dd6[k]]]] - 2*sqrt3*phi11;//	-->|_____|
 }							//		      A
 if(mask[dd2[k]]==1 && mask[dd3[k]]==1 && mask[dd5[k]]==1 && mask[dd2[dd3[dd5[k]]]]==28)	
 {							//		---|  x  |---
 p[dd2[dd3[dd5[k]]]] = p[dd1[dd4[dd6[k]]]] - 2*sqrt3*phi11;//	   |_____|<--
 }							//		      A
 if(mask[dd2[k]]==1 && mask[dd4[k]]==1 && mask[dd5[k]]==1 && mask[dd2[dd4[dd5[k]]]]==28)	
 {							//		---|  o  |---
 p[dd2[dd4[dd5[k]]]] = p[dd1[dd3[dd6[k]]]] - 2*sqrt3*phi11;//	   |_____|<--
 }							//		      A	
 */

/*
 //  top post edges
 if(mask[dd1[k]]==1 && mask[dd6[k]]==1 && mask[ddk[k]]==28)
 {											//				__v__
 p[ddk[k]] = p[ddj[k]] - 2*sqrt2*phi11;	    //edge A	-->|     |
 }											//			___|     |___
 if(mask[dd2[k]]==1 && mask[dd6[k]]==1 && mask[ddl[k]]==28)
 {											//				__v__	
 p[ddl[k]] = p[ddi[k]] - 2*sqrt2*phi11;	    //edge B	   |     |<--
 }											//			___|     |___	
 if(mask[dd3[k]]==1 && mask[dd6[k]]==1 && mask[ddg[k]]==28)
 {											//				__v__
 p[ddg[k]] = p[ddf[k]] - 2*sqrt2*phi11;	    //edge C	   |     |
 }											//			___|  x  |___
 if(mask[dd4[k]]==1 && mask[dd6[k]]==1 && mask[ddh[k]]==28)
 {											//				__v__
 p[ddh[k]] = p[dde[k]] - 2*sqrt2*phi11;	    //edge D	   |     |
 }											//			___|  o  |___
 
 
 //  around post edges
 if(mask[dd1[k]]==1 && mask[dd3[k]]==1 && mask[dda[k]]==28)
 {											//				_____
 p[dda[k]] = p[ddd[k]] - 2*sqrt2*phi11;	//edge E	-->|     |
 }											//			___|  x  |___
 if(mask[dd2[k]]==1 && mask[dd3[k]]==1 && mask[ddb[k]]==28)
 {											//				_____
 p[ddb[k]] = p[ddc[k]] - 2*sqrt2*phi11;	//edge F	   |     |<--
 }											//			___|  x  |___
 if(mask[dd4[k]]==1 && mask[dd2[k]]==1 && mask[ddd[k]]==28)
 {											//				_____
 p[ddd[k]] = p[dda[k]] - 2*sqrt2*phi11;	//edge G	   |     |<--	
 }											//			___|  o  |___
 if(mask[dd1[k]]==1 && mask[dd4[k]]==1 && mask[ddc[k]]==28)
 {											//				_____
 p[ddc[k]] = p[ddb[k]] - 2*sqrt2*phi11;	//edge H	-->|     |
 }											//			___|  o  |___
 
 
 
 // below posts edges
 if(mask[dd1[k]]==1 && mask[dd5[k]]==1 && mask[ddi[k]]==28)
 {											//			---|     |---
 p[ddi[k]] = p[ddl[k]] - 2*sqrt2*phi11;	//edge I	-->|_____|
 }											//		      A	
 if(mask[dd2[k]]==1 && mask[dd5[k]]==1 && mask[ddj[k]]==28)
 {											//			---|     |---
 p[ddj[k]] = p[ddk[k]] - 2*sqrt2*phi11;	//edge M	   |_____|<--
 }											//		      A
 if(mask[dd4[k]]==1 && mask[dd5[k]]==1 && mask[ddf[k]]==28)
 {											//			---|  o  |---
 p[ddf[k]] = p[ddg[k]] - 2*sqrt2*phi11;	//edge N	   |_____|
 }											//		      A
 if(mask[dd3[k]]==1 && mask[dd5[k]]==1 && mask[dde[k]]==28)
 {											//			---|  x  |---
 p[dde[k]] = p[ddh[k]] - 2*sqrt2*phi11;	//edge N	   |_____|
 }	 
 */





/*
// "ROUND" CORNERS
//  top post edges
if(mask[dd1[k]]==1 && mask[dd6[k]]==1 && mask[ddk[k]]==28)
{										
	p[ddk[k]] = factor*(p[ddi[k]] - 2*phi11 + p[ddl[k]] - 2*phi11) + factor2*(p[k] + 1/sqrt2*(p[ddj[k]]-p[k]) - (1+sqrt2)*phi11);
}	
if(mask[dd2[k]]==1 && mask[dd6[k]]==1 && mask[ddl[k]]==28)
{											
	p[ddl[k]] = factor*(p[ddj[k]] - 2*phi11 + p[ddk[k]] - 2*phi11) + factor2*(p[k] + 1/sqrt2*(p[ddi[k]]-p[k]) - (1+sqrt2)*phi11);			
}
if(mask[dd3[k]]==1 && mask[dd6[k]]==1 && mask[ddg[k]]==28)
{											
	p[ddg[k]] = factor*(p[dde[k]] - 2*phi11 + p[ddh[k]] - 2*phi11) + factor2*(p[k] + 1/sqrt2*(p[ddf[k]]-p[k]) - (1+sqrt2)*phi11);	
}
if(mask[dd4[k]]==1 && mask[dd6[k]]==1 && mask[ddh[k]]==28)
{	
	p[ddh[k]] = factor*(p[ddf[k]] - 2*phi11 + p[ddg[k]] - 2*phi11) + factor2*(p[k] + 1/sqrt2*(p[dde[k]]-p[k]) - (1+sqrt2)*phi11);				
}
//  around post edges
if(mask[dd1[k]]==1 && mask[dd3[k]]==1 && mask[dda[k]]==28)
{
	p[dda[k]] = factor*(p[ddb[k]] - 2*phi11 + p[ddc[k]] - 2*phi11) + factor2*(p[k] + 1/sqrt2*(p[ddd[k]]-p[k]) - (1+sqrt2)*phi11);	
}
if(mask[dd2[k]]==1 && mask[dd3[k]]==1 && mask[ddb[k]]==28)
{
	p[ddb[k]] = factor*(p[ddd[k]] - 2*phi11 + p[dda[k]] - 2*phi11) + factor2*(p[k] + 1/sqrt2*(p[ddc[k]]-p[k]) - (1+sqrt2)*phi11);	
}	
if(mask[dd4[k]]==1 && mask[dd2[k]]==1 && mask[ddd[k]]==28)
{	
	p[ddd[k]] = factor*(p[ddb[k]] - 2*phi11 + p[ddc[k]] - 2*phi11) + factor2*(p[k] + 1/sqrt2*(p[dda[k]]-p[k]) - (1+sqrt2)*phi11);	
}
if(mask[dd1[k]]==1 && mask[dd4[k]]==1 && mask[ddc[k]]==28)
{	
	p[ddc[k]] = factor*(p[ddd[k]] - 2*phi11 + p[dda[k]] - 2*phi11) + factor2*(p[k] + 1/sqrt2*(p[ddb[k]]-p[k]) - (1+sqrt2)*phi11);	
}

// below posts edges
if(mask[dd1[k]]==1 && mask[dd5[k]]==1 && mask[ddi[k]]==28)
{	
	p[ddi[k]] = factor*(p[ddj[k]] - 2*phi11 + p[ddk[k]] - 2*phi11) + factor2*(p[k] + 1/sqrt2*(p[ddl[k]]-p[k]) - (1+sqrt2)*phi11);	
}
if(mask[dd2[k]]==1 && mask[dd5[k]]==1 && mask[ddj[k]]==28)
{	
	p[ddj[k]] = factor*(p[ddi[k]] - 2*phi11 + p[ddl[k]] - 2*phi11) + factor2*(p[k] + 1/sqrt2*(p[ddk[k]]-p[k]) - (1+sqrt2)*phi11);	
}
if(mask[dd4[k]]==1 && mask[dd5[k]]==1 && mask[ddf[k]]==28)
{			      
	p[ddf[k]] = factor*(p[dde[k]] - 2*phi11 + p[ddh[k]] - 2*phi11) + factor2*(p[k] + 1/sqrt2*(p[ddg[k]]-p[k]) - (1+sqrt2)*phi11);	
}
if(mask[dd3[k]]==1 && mask[dd5[k]]==1 && mask[dde[k]]==28)
{	
	p[dde[k]] = factor*(p[ddf[k]] - 2*phi11 + p[ddg[k]] - 2*phi11) + factor2*(p[k] + 1/sqrt2*(p[ddh[k]]-p[k]) - (1+sqrt2)*phi11);		
}
*/
