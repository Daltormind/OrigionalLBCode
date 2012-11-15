//leakSearch.cpp: this function checks that all fluid near the surface is been counted
//Marco, 29 June 2011

#include "wet.h"
#include <iostream>

using namespace std;

void wet::leakSearch() 
{
	if (mask[k] == 0 && mask[dd1[k]] == 28) cout <<"leak"<< xk << " "<< yk <<" "<< zk << " " << 1 <<endl;
	if (mask[k] == 0 && mask[dd2[k]] == 28) cout <<"leak"<< xk << " "<< yk <<" "<< zk << " " << 2 <<endl;
	if (mask[k] == 0 && mask[dd3[k]] == 28) cout <<"leak"<< xk << " "<< yk <<" "<< zk << " " << 3 <<endl;
	if (mask[k] == 0 && mask[dd4[k]] == 28) cout <<"leak"<< xk << " "<< yk <<" "<< zk << " " << 4 <<endl;
	if (mask[k] == 0 && mask[dd5[k]] == 28) cout <<"leak"<< xk << " "<< yk <<" "<< zk << " " << 5 <<endl;
	if (mask[k] == 0 && mask[dd6[k]] == 28) cout <<"leak"<< xk << " "<< yk <<" "<< zk << " " << 6 <<endl;
	if (mask[k] == 0 && mask[dda[k]] == 28) cout <<"leak"<< xk << " "<< yk <<" "<< zk << " " << 7 <<endl;
	if (mask[k] == 0 && mask[ddb[k]] == 28) cout <<"leak"<< xk << " "<< yk <<" "<< zk << " " << 8 <<endl;
	if (mask[k] == 0 && mask[ddc[k]] == 28) cout <<"leak"<< xk << " "<< yk <<" "<< zk << " " << 9 <<endl;
	if (mask[k] == 0 && mask[ddd[k]] == 28) cout <<"leak"<< xk << " "<< yk <<" "<< zk << " " << 10 <<endl;
	if (mask[k] == 0 && mask[dde[k]] == 28) cout <<"leak"<< xk << " "<< yk <<" "<< zk << " " << 11 <<endl;
	if (mask[k] == 0 && mask[ddf[k]] == 28) cout <<"leak"<< xk << " "<< yk <<" "<< zk << " " << 12 <<endl;
	if (mask[k] == 0 && mask[ddg[k]] == 28) cout <<"leak"<< xk << " "<< yk <<" "<< zk << " " << 13 <<endl;
	if (mask[k] == 0 && mask[ddh[k]] == 28) cout <<"leak"<< xk << " "<< yk <<" "<< zk << " " << 14 <<endl;
	if (mask[k] == 0 && mask[ddi[k]] == 28) cout <<"leak"<< xk << " "<< yk <<" "<< zk << " " << 15 <<endl;
	if (mask[k] == 0 && mask[ddj[k]] == 28) cout <<"leak"<< xk << " "<< yk <<" "<< zk << " " << 16 <<endl;
	if (mask[k] == 0 && mask[ddk[k]] == 28) cout <<"leak"<< xk << " "<< yk <<" "<< zk << " " << 17 <<endl;
	if (mask[k] == 0 && mask[ddl[k]] == 28) cout <<"leak"<< xk << " "<< yk <<" "<< zk << " " << 18 <<endl;
}
