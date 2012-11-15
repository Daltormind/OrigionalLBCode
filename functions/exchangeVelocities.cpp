#include "wet.h"

void wet::exchangeVelocities(void)
{

	//if(t%infoStep==0)
	//	cout << "Process " << rank << ": exchanging phi density..." << endl;
	
	MPI_Status statusUxLeft, statusUxRight;
	MPI_Status statusUyLeft, statusUyRight;
	MPI_Status statusUzLeft, statusUzRight;
	
	MPI_Request requestOutUxLeft, requestOutUxRight, requestInUxLeft, requestInUxRight;
	MPI_Request requestOutUyLeft, requestOutUyRight, requestInUyLeft, requestInUyRight;
	MPI_Request requestOutUzLeft, requestOutUzRight, requestInUzLeft, requestInUzRight;
	

	MPI_Isend(&(uxs[k1]),   k1, MPI_DOUBLE, leftProcess,  rank*100 ,   MPI_COMM_WORLD, &requestOutUxLeft);
	MPI_Isend(&(uxs[k2-k1]),k1, MPI_DOUBLE, rightProcess, rank*100+1 , MPI_COMM_WORLD, &requestOutUxRight);
	MPI_Isend(&(uys[k1]),   k1, MPI_DOUBLE, leftProcess,  rank*100 ,   MPI_COMM_WORLD, &requestOutUyLeft);
	MPI_Isend(&(uys[k2-k1]),k1, MPI_DOUBLE, rightProcess, rank*100+1 , MPI_COMM_WORLD, &requestOutUyRight);
	MPI_Isend(&(uzs[k1]),   k1, MPI_DOUBLE, leftProcess,  rank*100 ,   MPI_COMM_WORLD, &requestOutUzLeft);
	MPI_Isend(&(uzs[k2-k1]),k1, MPI_DOUBLE, rightProcess, rank*100+1 , MPI_COMM_WORLD, &requestOutUzRight);
		
	MPI_Irecv(&(uxs[k2]), k1, MPI_DOUBLE, rightProcess, rightProcess*100, MPI_COMM_WORLD, &requestInUxLeft);
	MPI_Irecv(uxs,        k1, MPI_DOUBLE, leftProcess, leftProcess*100+1, MPI_COMM_WORLD, &requestInUxRight);
	MPI_Irecv(&(uys[k2]), k1, MPI_DOUBLE, rightProcess, rightProcess*100, MPI_COMM_WORLD, &requestInUyLeft);
	MPI_Irecv(uys,        k1, MPI_DOUBLE, leftProcess, leftProcess*100+1, MPI_COMM_WORLD, &requestInUyRight);
	MPI_Irecv(&(uzs[k2]), k1, MPI_DOUBLE, rightProcess, rightProcess*100, MPI_COMM_WORLD, &requestInUzLeft);
	MPI_Irecv(uzs,        k1, MPI_DOUBLE, leftProcess, leftProcess*100+1, MPI_COMM_WORLD, &requestInUzRight);
	
	//Waiting until ux are saved in the recieving buffer
	MPI_Wait(&requestInUxLeft,  &statusUxLeft);
	MPI_Wait(&requestInUxRight, &statusUxRight);
	MPI_Wait(&requestInUyLeft,  &statusUyLeft);
	MPI_Wait(&requestInUyRight, &statusUyRight);
	MPI_Wait(&requestInUzLeft,  &statusUzLeft);
	MPI_Wait(&requestInUzRight, &statusUzRight);

	
	//Waiting until sending buffer is relased
	MPI_Wait(&requestOutUxLeft, &statusUxLeft);
	MPI_Wait(&requestOutUxRight,&statusUxRight);
	MPI_Wait(&requestOutUyLeft, &statusUyLeft);
	MPI_Wait(&requestOutUyRight,&statusUyRight);
	MPI_Wait(&requestOutUzLeft, &statusUzLeft);
	MPI_Wait(&requestOutUzRight,&statusUzRight);

		

	
	//if(t%infoStep==0)
	//	cout << "Process " << rank << ": phi density exchanged." << endl;	
}
