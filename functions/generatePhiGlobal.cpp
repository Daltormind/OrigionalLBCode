#include "wet.h"

void wet::generatePhiGlobal(void)
{

	//if(t%infoStep==0)
	//	cout << "Process " << rank << ": generating global phi..." << endl;

	int n, m, g;
	MPI_Status  status[size];
	
	MPI_Request request[size];


	if(t%writeStep==0)
	{	
		if(rank!=ROOT)
			MPI_Isend(&(p[k1]), k2-k1, MPI_DOUBLE, ROOT, rank, MPI_COMM_WORLD, &request[rank]);
		if(rank==ROOT)
		{
			for(m=0; m<size; m++)
			{
				if(m<ROOT)
					MPI_Irecv(&(pGlobal[m*(k2-k1*(1+LX%size))]), k2-k1*(1+LX%size), MPI_DOUBLE, m, m, MPI_COMM_WORLD, &request[m]);
				if(m==ROOT)
				{
					for(g=0; g<k2-k1; g++)
						pGlobal[m*(k2-k1*(1+LX%size))+g]=p[g+k1];
				}
				if(m>ROOT)
					MPI_Irecv(&(pGlobal[m*(k2-k1*(1+LX%size))+(LX%size)*k1]), k2-k1*(1+LX%size), MPI_DOUBLE, m, m, 	MPI_COMM_WORLD, &request[m]);
			}
		}
		
		
		if(rank==ROOT)				
		{
			for(n=0; n<size; n++)
			{
				if(n!=ROOT)
					MPI_Wait(&request[n], &status[n]);
			}
		}
	
	}

	//if(t%infoStep==0)
	//	cout << "Process " << rank << ": global phi generated." << endl;
}
			
