//writeYPlanXVelocityFile(Ylayer).cpp: for a certain Y-layer it writes the x and z velocity components into matlab m-file: uxt(:,:)= [...]
//filename: uxYytT.m where T = number of iterations and y is the Y-layer number
//Lisa, 28 July 2011

#include "wet.h"



void wet::writeDissipationFile(void)
{
	
	int i,j,h;
	 
	
	//cout << " Process " << rank << ": check1 " << endl;
	//dissGlobal = generateGlobal(dissipation);
	//cout << " Process " << rank << ": check2 " << endl;
	
	
	
	
	MPI_Gather(&(dissipation[k1]),k2-k1,MPI_DOUBLE,dissGlobal    ,k2-k1,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	MPI_Gather(&(diffDiss[k1])   ,k2-k1,MPI_DOUBLE,diffDissGlobal,k2-k1,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	//MPI_Gather(&(eBalDiff[k1])   ,k2-k1,MPI_DOUBLE,eBalDiffGlobal,k2-k1,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);

	if(rank==ROOT)
	{
		int i, j, h;	
		
		char filename[25];
		sprintf(filename, "./densitydata/ddiss%ld.m", t);			//Create a name for file that contain data 
		//sprintf(filename, "./dd%ld.m", t);
		
		ofstream file(filename);
		file.precision(4);
		for( h = 0 ; h < LZ ; h++) 
		{   
			file << "diss" << t << "(:,:," << h+1 << ")=[" << endl;
			for( i = 0 ; i < LX ; i++) 
			{
				for( j = 0 ; j < LY ; j++) 
				{
					k = h + j*LZ + i*LY*LZ;
					if( maskGlobal[k]==28) 
						file << 0.0 << " ";  
					else
						file << dissGlobal[k] << " ";
					
				}
				file << endl;
			}
			file <<"];" << endl;
		}
		file.close();
		
		char filename2[25];
		sprintf(filename2, "./densitydata/ddiffdiss%ld.m", t);			//Create a name for file that contain data 

		ofstream file2(filename2);
		file2.precision(4);
		for( h = 0 ; h < LZ ; h++) 
		{   
			file2 << "diffdiss" << t << "(:,:," << h+1 << ")=[" << endl;
			for( i = 0 ; i < LX ; i++) 
			{
				for( j = 0 ; j < LY ; j++) 
				{
					k = h + j*LZ + i*LY*LZ;
					if( maskGlobal[k]==28) 
						file2 << 0.0 << " ";  
					else
						file2 << diffDissGlobal[k] << " ";
					
				}
				file2 << endl;
			}
			file2 <<"];" << endl;
		}
		file2.close();
		
		
		/*
		char filename3[25];
		sprintf(filename3, "./densitydata/deBal%ld.m", t);			//Create a name for file that contain data 
		
		ofstream file3(filename3);
		file3.precision(4);
		for( h = 0 ; h < LZ ; h++) 
		{   
			file3 << "eBal" << t << "(:,:," << h+1 << ")=[" << endl;
			for( i = 0 ; i < LX ; i++) 
			{
				for( j = 0 ; j < LY ; j++) 
				{
					k = h + j*LZ + i*LY*LZ;
					if( maskGlobal[k]==28) 
						file3 << 0.0 << " ";  
					else
						file3 << eBalDiffGlobal[k] << " ";
					
				}
				file3 << endl;
			}
			file3 <<"];" << endl;
		}
		file3.close();
		 */
	}

	
}

