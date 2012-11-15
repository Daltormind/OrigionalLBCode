//saveFiles.cpp: save data in files
//Marco, 30 June 2011, Lisa, 28 July 2011

#include "wet.h"



void wet::saveFiles(void)
{	
		
	
	if (rank == ROOT){
		if(t%(infoStep) == 0){
			{
                char filenameE[50];
                snprintf(filenameE,50, "./uxs%.3frat%.3foff%.0f/energy.dat",initUX,dropletR/Dx,dropletCenterY-Px-Dx/2);
				ofstream file20;
				file20.open(filenameE, ios::app);
				file20.precision(12);
				file20 << t << " " << energy << " " << bulkE<< " " << interfaceE << " " << surfaceE << "      " << energy_n << " " << bulkE_n << " " << interfaceE_n << " " << surfaceE_n << "      " << kin << " " << kin_CM <<  endl;   
				file20.close();
				
			}
		}
	}
	if(t%writeStep == 0 )
	{	
		writeDensityFile();
		
		if (rank == ROOT) cout << "time = "<< t << ", total free energy " << energy_n << " , contact Area " << surfArea <<  endl;
			
		
		if (velocityprofile == true) {
			//writeYPlanXVelocityFile(Dy/2);
			writeZPlanXVelocityFile(1);
			//writeYPlanXVelocityFile( (LY+Dy)/2);
		}
		//cout << "Process " << rank << ": saveFiles check 3, t = "<< t << endl;
		//writeDissipationFile(); Removed as not currently interested in disipation Mat 18/10/2012
		//cout << "Process " << rank << ": saveFiles check 4, t = "<< t << endl;
		
	}
	
	/*if(t % (infoStep) == 0)
	{
		energy = bulkE + interfaceE + surfaceE;
		ofstream file;
		//file.open("./data/energy.dat", ios::app);
		file.open("./energy.dat", ios::app);
		file.precision(12);
		file << "t= " << t << " energy= " << energy << " bulk-energy= " << bulkE<< " interface-energy= " << interfaceE << "surface-energy= " << surfaceE << endl;
		file.close();
		
	}*/
	
}

