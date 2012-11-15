#include "wet.h"

void wet::exchangeMask(void)
{

	//if(MPI_Barrier(MPI_COMM_WORLD)==MPI_SUCCESS)
	//	cout << "Process " << rank << ": exchanging mask value between processes ring..." << endl;

	MPI_Status statusLeft, statusRight;
	
	//Left part sendind and right part reading (using blocking function=slow but safe (THINK SOMETHING BETTER))
	if(rank==ROOT)
	{
		MPI_Send(&(mask[k1]), k1, MPI_INT, leftProcess, rank, MPI_COMM_WORLD);
		MPI_Recv(&(mask[k2]), k1, MPI_INT, rightProcess, rightProcess, MPI_COMM_WORLD, &statusRight); 
	}
	else
	{
		MPI_Recv(&(mask[k2]), k1, MPI_INT, rightProcess, rightProcess, MPI_COMM_WORLD, &statusLeft); 
		MPI_Send(&(mask[k1]), k1, MPI_INT, leftProcess, rank, MPI_COMM_WORLD);
	}
	
	//Right part sending and left part reciving (using blocking function=slow but safe (THINK SOMETHING BETTER))
	if(rank==ROOT)
	{	
		MPI_Send(&(mask[k2-k1]), k1, MPI_INT, rightProcess, rank, MPI_COMM_WORLD);
		MPI_Recv(mask, k1, MPI_INT, leftProcess, leftProcess, MPI_COMM_WORLD, &statusLeft);
	}
	
	else
	{
		MPI_Recv(mask, k1, MPI_INT, leftProcess, leftProcess, MPI_COMM_WORLD, &statusLeft);
		MPI_Send(&(mask[k2-k1]), k1, MPI_INT, rightProcess, rank, MPI_COMM_WORLD);
	}
	//	cout << "Process " << rank << ": mask values exchanged." << endl;
}
