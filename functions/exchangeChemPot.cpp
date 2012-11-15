#include "wet.h"

void wet::exchangeChemPot(void)
{

	//if(t%infoStep==0)
	//	cout << "Process " << rank << ": exchanging phi density..." << endl;
	
	MPI_Status statusMuLeft, statusMuRight;
	
	MPI_Request requestOutMuLeft, requestOutMuRight, requestInMuLeft, requestInMuRight;
	
	MPI_Isend(&(mu[k1]),   k1, MPI_DOUBLE, leftProcess,  rank*100 ,   MPI_COMM_WORLD, &requestOutMuLeft);
	MPI_Isend(&(mu[k2-k1]),k1, MPI_DOUBLE, rightProcess, rank*100+1 , MPI_COMM_WORLD, &requestOutMuRight);
		
	MPI_Irecv(&(mu[k2]), k1, MPI_DOUBLE, rightProcess, rightProcess*100, MPI_COMM_WORLD, &requestInMuLeft);
	MPI_Irecv(mu,        k1, MPI_DOUBLE, leftProcess, leftProcess*100+1, MPI_COMM_WORLD, &requestInMuRight);
	
	//Waiting until mu are saved in the recieving buffer
	MPI_Wait(&requestInMuLeft,  &statusMuLeft);
	MPI_Wait(&requestInMuRight, &statusMuRight);
		
	//Waiting until sending buffer is relased
	MPI_Wait(&requestOutMuLeft, &statusMuLeft);
	MPI_Wait(&requestOutMuRight,&statusMuRight);
	

}
