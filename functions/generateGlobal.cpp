// generateGlobal.cpp: Takes an array (sub-arrays of size processN from all CPUs) 
// and combines them in a global array of size N on the ROOT CPU
// Lisa, 02 Aug 2011

#include "wet.h"

double* wet::generateGlobal(double *array)
{

	//if(t%infoStep==0)
	//	cout << "Process " << rank << ": generating global n..." << endl;

	int m, g;
	double *arrayGlobal;
	
	arrayGlobal = new double[N]; 
	
	MPI_Status  status[size];
	MPI_Request request[size];



		
	if(rank!=ROOT)
		MPI_Isend(&array[k1], k2-k1, MPI_DOUBLE, ROOT, rank, MPI_COMM_WORLD, &request[rank]);
	if(rank==ROOT)
	{
		for(m=0; m<size; m++)
		{
			if(m<ROOT)
				MPI_Irecv(&(arrayGlobal[m*(k2-k1*(1+LX%size))]), k2-k1*(1+LX%size), MPI_DOUBLE, m, m, MPI_COMM_WORLD, &request[m]);
			if(m==ROOT)
			{
				for(g=0; g<k2-k1; g++)
					arrayGlobal[m*(k2-k1*(1+LX%size))+g]=array[g+k1];
			}
			if(m>ROOT)
				MPI_Irecv(&(arrayGlobal[m*(k2-k1*(1+LX%size))+(LX%size)*k1]), k2-k1*(1+LX%size), MPI_DOUBLE, m, m, 	MPI_COMM_WORLD, &request[m]);
		}
	}
	
	
	if(rank==ROOT)				
	{
		for(m=0; m<size; m++)
		{
			if(m!=ROOT)
				MPI_Wait(&request[m], &status[m]);
		}
	}
	
	
	
	return arrayGlobal;
	//if(t%infoStep==0)
	//	cout << "Process " << rank << ": global n generated." << endl;
}
			
