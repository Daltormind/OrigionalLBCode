//calculatePenetrationLength.cpp: 
//Marco, 30 June 2011

#include "wet.h"
#include <fstream>
using namespace std;


void wet::calculatePenetrationLength(void)
{

	ofstream file("./data/Length.m");
	file.precision(4);
	
	j = 0; h = LZ/2;
	for ( i = 0; i < length; i++ ) 
	{
		if (p[h + j*LZ + i*LY*LZ] * p[h + j*LZ + (i+1)*LY*LZ] <= 0 && p[h + j*LZ + (i+1)*LY*LZ] != 0.0) 
		{
			file << t << i + ( p[h + j*LZ + i*LY*LZ] ) / ( p[h + j*LZ + i*LY*LZ] - p[h + j*LZ + (i+1)*LY*LZ] ) << endl;
		}
	}
	file.close();
}
