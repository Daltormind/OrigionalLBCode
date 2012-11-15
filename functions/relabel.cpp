//relabel.cpp: this function rename the mask of fluid (mask=0) near a solid surface with mask=1, and the solidDist from 0 to 0.5
//Marco, 29 June 2011

#include "wet.h"

void wet::relabel(void)
{
	
	if (mask[k] == 0) 
	{

		if (mask[dd1[k]] == 28)
			mask[k] = 1;
		if (mask[dd2[k]] == 28)
			mask[k] = 1;
		if (mask[dd3[k]] == 28)
			mask[k] = 1;
		if (mask[dd4[k]] == 28)
			mask[k] = 1;
		if (mask[dd5[k]] == 28) 
			mask[k] = 1;
		if (mask[dd6[k]] == 28)
			mask[k] = 1;
		if (mask[dda[k]] == 28)
			mask[k] = 1;
		if (mask[ddb[k]] == 28)
			mask[k] = 1;
		if (mask[ddc[k]] == 28)
			mask[k] = 1;
		if (mask[ddd[k]] == 28)
			mask[k] = 1;
		if (mask[dde[k]] == 28)
			mask[k] = 1;
		if (mask[ddf[k]] == 28)
			mask[k] = 1;
		if (mask[ddg[k]] == 28)
			mask[k] = 1;
		if (mask[ddh[k]] == 28)
			mask[k] = 1;
		if (mask[ddi[k]] == 28)
			mask[k] = 1;
		if (mask[ddj[k]] == 28)
			mask[k] = 1;
		if (mask[ddk[k]] == 28)
			mask[k] = 1;
		if (mask[ddl[k]] == 28)
			mask[k] = 1;
	}
	
}

