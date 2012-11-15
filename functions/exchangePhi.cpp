#include "wet.h"

void wet::exchangePhi(void)
{

	//if(t%infoStep==0)
	//	cout << "Process " << rank << ": exchanging phi density..." << endl;
	
	MPI_Status statusPhiLeft, statusPhiRight;
	
	MPI_Request requestOutPhiLeft, requestOutPhiRight, requestInPhiLeft, requestInPhiRight;
	
	//Sending left (for sender) phi total density
	MPI_Isend(&(p[k1]),k1, MPI_DOUBLE, leftProcess, rank*100 , MPI_COMM_WORLD, &requestOutPhiLeft);
	//Sending right (for sender) phi total density
	MPI_Isend(&(p[k2-k1]),k1, MPI_DOUBLE, rightProcess, rank*100+1 , MPI_COMM_WORLD, &requestOutPhiRight);
		
	//Recieving right (for reciever) phi total density
	MPI_Irecv(&(p[k2]), k1, MPI_DOUBLE, rightProcess, rightProcess*100, MPI_COMM_WORLD, &requestInPhiLeft);
	//Recieving left (for reciever) phi total density
	MPI_Irecv(p, k1, MPI_DOUBLE, leftProcess, leftProcess*100+1, MPI_COMM_WORLD, &requestInPhiRight);
	
	//Waiting until phi densities are saved in the recieving buffer
	MPI_Wait(&requestInPhiLeft, &statusPhiLeft);
	MPI_Wait(&requestInPhiRight, &statusPhiRight);
	
	//Waiting until sending buffer is relased
	MPI_Wait(&requestOutPhiLeft, &statusPhiLeft);
	MPI_Wait(&requestOutPhiRight, &statusPhiRight);
	
	
	
	//Sending left (for sender) n total density
	MPI_Isend(&(n[k1]),k1, MPI_DOUBLE, leftProcess, rank*100 , MPI_COMM_WORLD, &requestOutPhiLeft);
	//Sending right (for sender) phi total density
	MPI_Isend(&(n[k2-k1]),k1, MPI_DOUBLE, rightProcess, rank*100+1 , MPI_COMM_WORLD, &requestOutPhiRight);
	
	//Recieving right (for reciever) n total density
	MPI_Irecv(&(n[k2]), k1, MPI_DOUBLE, rightProcess, rightProcess*100, MPI_COMM_WORLD, &requestInPhiLeft);
	//Recieving left (for reciever) phi total density
	MPI_Irecv(n, k1, MPI_DOUBLE, leftProcess, leftProcess*100+1, MPI_COMM_WORLD, &requestInPhiRight);
	
	//Waiting until n densities are saved in the recieving buffer
	MPI_Wait(&requestInPhiLeft, &statusPhiLeft);
	MPI_Wait(&requestInPhiRight, &statusPhiRight);
	
	//Waiting until sending buffer is relased
	MPI_Wait(&requestOutPhiLeft, &statusPhiLeft);
	MPI_Wait(&requestOutPhiRight, &statusPhiRight);
		

	
	//if(t%infoStep==0)
	//	cout << "Process " << rank << ": phi density exchanged." << endl;	
}
