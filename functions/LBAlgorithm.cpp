//LBAlgorithm.cpp: lattice boltzmann algorithm
//Marco, 29 June 2011, Lisa, 28 July 2011

#include "wet.h"


void wet::LBAlgorithm(int startStep, int runStep)
{	
	
	cout << "- in LBAlgorithm ------------ " <<endl;
	cout << "runStep " << runStep << endl;
	cout << "nbEqStep " << nbEqStep << endl;
	cout << "equilTime " << equilTime << endl;
	cout << "teta1 " << teta1 << endl;
	cout << "tauliquid " << tauliquid << endl;
	cout << "----------------------------- " <<endl;
	
	
	//cout << "Process " << rank << ": LBAlgorithm check 0" << endl;
	
    //for(t = startStep; t <= runStep; t++)
    do
	{		
		
		//cout << "Process " << rank << ": , t = " << t << " LBAlgorithm check 1" << endl;
		if (t%infoStep==0)
		{
			//variables needed for computeFreeEnergy
			energy = 0.0;
			bulkE = 0.0;
			interfaceE = 0.0;
			surfaceE = 0.0;
			energy_n = 0.0;
			bulkE_n = 0.0;
			interfaceE_n = 0.0;
			surfaceE_n = 0.0;
		
		}
        
		//cout << "Process " << rank << ": , t = " << t << " LBAlgorithm check 1, " << "uzs[k1] "<< uxs[k1] << endl;
		
		//cout << "Process " << rank << ": , t = " << t << " LBAlgorithm check 1, ff0[k1+Dh] "<< ff0[k1+Dh] << " , ff0[k1+Dh+LY*LZ] "<< ff0[k1+Dh+LY*LZ] << " , ff0[k1+Dh+2*LY*LZ] "<< ff0[k1+Dh+2*LY*LZ]<< endl;		
		//cout << "Process " << rank << ": , t = " << t << " LBAlgorithm check 2, " << "uzs[k1] "<< uxs[k1] << endl;
		/*exchangePhi();//**
		exchangeDensities();//**
		exchangeDensities_ffgg();//**
		exchangeVelocities();//**
		*/
		
		//cout << "Process " << rank << ": , t = " << t << " LBAlgorithm check 2, ff0[k1+Dh] "<< ff0[k1+Dh] << " , ff0[k1+Dh+LY*LZ] "<< ff0[k1+Dh+LY*LZ] << " , ff0[k1+Dh+2*LY*LZ] "<< ff0[k1+Dh+2*LY*LZ]<< endl;

		//cout << "Process " << rank << ": , t = " << t << " LBAlgorithm check 2, uzs[k1+Dh] "<< uzs[k1+Dh] << " , uzs[k1+Dh+LY*LZ] "<< uzs[k1+Dh+LY*LZ] << endl;
		
        
        computeMomenta();
		
        //cout << "Process " << rank << ": , t = " << t << " LBAlgorithm check 3, uzs[k1+Dh] "<< uzs[k1+Dh] << " , uzs[k1+Dh+LY*LZ] "<< uzs[k1+Dh+LY*LZ] << endl;
		
		
		
		
		exchangePhi();			//here to calculate densities at solid boundaries
		
				
		//cout << "Process " << rank << ": , t = " << t << " LBAlgorithm check 4, " << "p[k1] "<< p[k1] << endl;
		densitiesAtSolidBoundaries();

		//cout << "Process " << rank << ": , t = " << t << " LBAlgorithm check 5, " << "p[k1] "<< p[k1] << endl;
				
		exchangePhi();			//and here to calculate equilibrium() function in collition()
		//cout << "Process " << rank << ": , t = " << t << " LBAlgorithm check 6, " << "p[k1] "<< p[k1] << endl;
		
		
		
		collision(); 
		//cout << "Process " << rank << ": , t = " << t << " LBAlgorithm check 7, " << "uxs[k1] "<< uxs[k1] << endl;
		
		
		
		generatePhiGlobal();
		//cout << "Process " << rank << ": , t = " << t << " LBAlgorithm check 8, " << "p[k1] "<< p[k1] << endl;
		
		cout << "Got to just before velav rank =" << rank << endl;
		
		
		velav();
		
		cout << "Got to just agter velav rank =" << rank << endl;
		if(t%infoStep==0) {
			ARolling();
			computeContactArea();
			exchangeVelocities();
			exchangeChemPot();
			computeDissipation();
            			
		}
		
		saveFiles();
		
		
				
		//cout << "Process " << rank << ": , t = " << t << " LBAlgorithm check 10, " << "p[k1] "<< p[k1] << endl;
		exchangeDensities();
		
				//cout << "Process " << rank << ": , t = " << t << " LBAlgorithm check 10, " << "uxs[k1] "<< uxs[k1] << endl;
		propagation();
		/*
		exchangePhi();//**
		exchangeDensities();//**
		exchangeDensities_ffgg();//**
		exchangeVelocities();//**
		 */
		//cout << "Process " << rank << ": , t = " << t << " LBAlgorithm check 11, " << "uxs[k1] "<< uxs[k1] << endl;
		exchangeDensities_ffgg();
		applyBoundaryConditions();
		exchangeDensities_ffgg();
		/*
		exchangePhi();//**
		exchangeDensities();//**
		exchangeDensities_ffgg();//**
		exchangeVelocities();//**	
		 */
		//cout << "Process " << rank << ": , t = " << t << " LBAlgorithm check 12, " << "uxs[k1] "<< uxs[k1] << endl;
		//MPI_Barrier(MPI_COMM_WORLD);
	t+=1;
	
	cout << "co is erqual to" << co << endl;
	}while(co<10);
}

