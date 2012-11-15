//writeDensityFile.cpp: writes the concentration into matlab m-file: dt(:,:,h+1)= [...]
//filename: ddt.m where t = number of iterations and h is zk (i.e. density is written in x-y plane slices)
//Marco, 30 June 2011, Lisa, 28 July 2011

#include "wet.h"

void wet::writeDensityFile(void)
{  
	
	if(rank==ROOT)
	{
		int i, j, h;	

		char filename[50];
		sprintf(filename, "./uxs%.3frat%.3foff%.0f/dd%ld.m",initUX,dropletR/Dx,dropletCenterY-Px-Dx/2, t);			//Create a name for file that contain data
		//sprintf(filename, "./dd%ld.m", t);
		
		ofstream file(filename);
		file.precision(4);
		for( h = 0 ; h < LZ ; h++) 
		{   
			file << "d" << t << "(:,:," << h+1 << ")=[" << endl;
			for( i = 0 ; i < LX ; i++) 
			{
				for( j = 0 ; j < LY ; j++) 
				{
					k = h + j*LZ + i*LY*LZ;
					if( maskGlobal[k]==28) 
						file << -2.0 << " ";  
					else
						file << pGlobal[k] << " ";
						
				}
				file << endl;
			}
			file <<"];" << endl;
		}
		file.close();
	}
	
	/*
	//generateNGlobal();	//ACHTUNG! this one causes the calculation on NEWHYDRA to crash!
	if(rank==ROOT)
	{
		int i, j, h;	
		
		char filename2[25];
		sprintf(filename2, "./densitydata/dn%ld.m", t);			//Create a name for file that contain data 
		//sprintf(filename, "./dd%ld.m", t);
		
		ofstream file(filename2);
		file.precision(4);
		for( h = 0 ; h < LZ ; h++) 
		{   
			file << "n" << t << "(:,:," << h+1 << ")=[" << endl;
			for( i = 0 ; i < LX ; i++) 
			{
				for( j = 0 ; j < LY ; j++) 
				{
					k = h + j*LZ + i*LY*LZ;
					if( maskGlobal[k]==28) 
						file << -2.0 << " ";  
					else
						file << nGlobal[k] << " ";
					
				}
				file << endl;
			}
			file <<"];" << endl;
		}
		file.close();
	}
	*/

/*
	muGlobal=generateGlobal(mu);
	if(rank==ROOT)
	{
		int i, j, h;	
		
		char filename3[25];
		sprintf(filename3, "./densitydata/dmu%ld.m", t);			//Create a name for file that contain data 
		//sprintf(filename, "./dd%ld.m", t);
		
		ofstream file(filename3);
		file.precision(4);
		for( h = 0 ; h < LZ ; h++) 
		{   
			file << "mu" << t << "(:,:," << h+1 << ")=[" << endl;
			for( i = 0 ; i < LX ; i++) 
			{
				for( j = 0 ; j < LY ; j++) 
				{
					k = h + j*LZ + i*LY*LZ;
					if( maskGlobal[k]==28) 
						file << 0 << " ";  
					else
						file << muGlobal[k] << " ";
					
				}
				file << endl;
			}
			file <<"];" << endl;
		}
		file.close();
	}
	
	*/


}

