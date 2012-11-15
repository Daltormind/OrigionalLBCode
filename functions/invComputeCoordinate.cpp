//invComputeCoordinate.cpp: calculate the k index from the node indeces
//Marco, 30 June 2011

#include "wet.h"
using namespace std;

void wet::invComputeCoordinate(void)                     
{                                             
	k = xk*LZ*LY + yk*LZ + zk;
}
