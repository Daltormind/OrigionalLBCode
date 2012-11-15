#include "wet.h"

void wet::generateNGlobal(void)
{

	//if(t%infoStep==0)
	//	cout << "Process " << rank << ": generating global n..." << endl;

	int m, g;
	MPI_Status  status[size];
	
	MPI_Request request[size];


	if(t%writeStep==0)
	{	
		if(rank!=ROOT)
			MPI_Isend(&(n[k1]), k2-k1, MPI_DOUBLE, ROOT, rank, MPI_COMM_WORLD, &request[rank]);
		if(rank==ROOT)
		{
			for(m=0; m<size; m++)
			{
				if(m<ROOT)
					MPI_Irecv(&(nGlobal[m*(k2-k1*(1+LX%size))]), k2-k1*(1+LX%size), MPI_DOUBLE, m, m, MPI_COMM_WORLD, &request[m]);
				if(m==ROOT)
				{
					for(g=0; g<k2-k1; g++)
						nGlobal[m*(k2-k1*(1+LX%size))+g]=n[g+k1];
				}
				if(m>ROOT)
					MPI_Irecv(&(nGlobal[m*(k2-k1*(1+LX%size))+(LX%size)*k1]), k2-k1*(1+LX%size), MPI_DOUBLE, m, m, 	MPI_COMM_WORLD, &request[m]);
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
	
	}

	//if(t%infoStep==0)
	//	cout << "Process " << rank << ": global n generated." << endl;
}
			
