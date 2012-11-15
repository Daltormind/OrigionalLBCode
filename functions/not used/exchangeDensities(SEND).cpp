#include "wet.h"

void wet::exchangeDensities(void)
{
	if(t%infoStep==0)
		cout << "Process " << rank << ": exchanging densities...." << endl;
	
	MPI_Status statusfn1, statusfn2, statusfna, statusfnb, statusfnc, statusfnd, statusfni, statusfnj, statusfnk, statusfnl, statusgn1, statusgn2, statusgna, statusgnb, statusgnc, statusgnd, statusgni, statusgnj, statusgnk, statusgnl; 
	
	
	MPI_Request requestOutfn1, requestOutfn2, requestInfn1, requestInfn2, requestOutfna, requestOutfnb, requestInfna, requestInfnb, requestOutfnc, requestOutfnd, requestInfnc, requestInfnd, requestOutfni, requestOutfnj, requestInfni, requestInfnj, requestOutfnk, requestOutfnl, requestInfnk, requestInfnl,requestOutgn1, requestOutgn2, requestIngn1, requestIngn2, requestOutgna, requestOutgnb, requestIngna, requestIngnb, requestOutgnc, requestOutgnd, requestIngnc, requestIngnd, requestOutgni, requestOutgnj, requestIngni, requestIngnj, requestOutgnk, requestOutgnl, requestIngnk, requestIngnl;


	if(rank==ROOT)
	{
	
	//Sending right (for sender) fn1 (only right part, because fn1 moves in x direction)
	MPI_Send(&(fn1[k2-k1]),k1, MPI_DOUBLE, rightProcess, rank*2 , MPI_COMM_WORLD);
	//Recieving left (for reciever) fn1 (only left part, because fn1 moves in x direction)
	MPI_Recv(fn1, k1, MPI_DOUBLE, leftProcess, leftProcess*2, MPI_COMM_WORLD, &statusfn1);
	
	//Sending left (for sender) fn2 (only left part, because fn2 moves in -x direction)
	MPI_Send(&(fn2[k1]),k1, MPI_DOUBLE, leftProcess, rank*3 , MPI_COMM_WORLD);
	//Recieving right (for reciever) fn2 (only right part, because fn2 moves in -x direction)
	MPI_Recv(&(fn2[k2]), k1, MPI_DOUBLE, rightProcess, rightProcess*3, MPI_COMM_WORLD, &statusfn2);

	//Sending right (for sender) fna (only right part, because fna moves in x direction)
	MPI_Send(&(fna[k2-k1]),k1, MPI_DOUBLE, rightProcess, rank*5 , MPI_COMM_WORLD);
	//Recieving left (for reciever) fna (only left part, because fna moves in x direction)
	MPI_Recv(fna, k1, MPI_DOUBLE, leftProcess, leftProcess*5, MPI_COMM_WORLD, &statusfna);

	//Sending left (for sender) fnb (only left part, because fnb moves in -x direction)
	MPI_Send(&(fnb[k1]),k1, MPI_DOUBLE, leftProcess, rank*7 , MPI_COMM_WORLD);
	//Recieving right (for reciever) fnb (only right part, because fnb moves in -x direction)
	MPI_Recv(&(fn2[k2]), k1, MPI_DOUBLE, rightProcess, rightProcess*7, MPI_COMM_WORLD, &statusfnb);

	//Sending right (for sender) fnc (only right part, because fnc moves in x direction)
	MPI_Send(&(fnc[k2-k1]),k1, MPI_DOUBLE, rightProcess, rank*11 , MPI_COMM_WORLD);
	//Recieving left (for reciever) fn1 (only left part, because fnc moves in x direction)
	MPI_Recv(fnc, k1, MPI_DOUBLE, leftProcess, leftProcess*11, MPI_COMM_WORLD, &statusfnc);

	//Sending left (for sender) fnd (only left part, because fnd moves in -x direction)
	MPI_Send(&(fnd[k1]),k1, MPI_DOUBLE, leftProcess, rank*13 , MPI_COMM_WORLD);
	//Recieving right (for reciever) fn2 (only right part, because fnd moves in -x direction)
	MPI_Recv(&(fnd[k2]), k1, MPI_DOUBLE, rightProcess, rightProcess*13, MPI_COMM_WORLD, &statusfnd);

	//Sending right (for sender) fni (only right part, because fni moves in x direction)
	MPI_Send(&(fni[k2-k1]),k1, MPI_DOUBLE, rightProcess, rank*17 , MPI_COMM_WORLD);
	//Recieving left (for reciever) fna (only left part, because fni moves in x direction)
	MPI_Recv(fni, k1, MPI_DOUBLE, leftProcess, leftProcess*17, MPI_COMM_WORLD, &statusfni);
		
	//Sending left (for sender) fnj (only left part, because fnj moves in -x direction)
	MPI_Send(&(fnj[k1]),k1, MPI_DOUBLE, leftProcess, rank*19 , MPI_COMM_WORLD);
	//Recieving right (for reciever) fnb (only right part, because fnj moves in -x direction)
	MPI_Recv(&(fnj[k2]), k1, MPI_DOUBLE, rightProcess, rightProcess*19, MPI_COMM_WORLD, &statusfnj);

	//Sending right (for sender) fnk (only right part, because fnk moves in x direction)
	MPI_Send(&(fnk[k2-k1]),k1, MPI_DOUBLE, rightProcess, rank*23 , MPI_COMM_WORLD);
	//Recieving left (for reciever) fnk (only left part, because fnk moves in x direction)
	MPI_Recv(fnk, k1, MPI_DOUBLE, leftProcess, leftProcess*23, MPI_COMM_WORLD, &statusfnk);

	//Sending left (for sender) fnl (only left part, because fnl moves in -x direction)
	MPI_Send(&(fnl[k1]),k1, MPI_DOUBLE, leftProcess, rank*29 , MPI_COMM_WORLD);
	//Recieving right (for reciever) fnl (only right part, because fnl moves in -x direction)
	MPI_Recv(&(fnl[k2]), k1, MPI_DOUBLE, rightProcess, rightProcess*29, MPI_COMM_WORLD, &statusfnl);

	//Sending right (for sender) gn1 (only right part, because gn1 moves in x direction)
	MPI_Send(&(gn1[k2-k1]),k1, MPI_DOUBLE, rightProcess, rank*31 , MPI_COMM_WORLD);
	//Recieving left (for reciever) gn1 (only left part, because gn1 moves in x direction)
	MPI_Recv(gn1, k1, MPI_DOUBLE, leftProcess, leftProcess*31, MPI_COMM_WORLD, &statusgn1);

	//Sending left (for sender) gn2 (only left part, because gn2 moves in -x direction)
	MPI_Send(&(gn2[k1]),k1, MPI_DOUBLE, leftProcess, rank*37 , MPI_COMM_WORLD);
	//Recieving right (for reciever) gn2 (only right part, because gn2 moves in -x direction)
	MPI_Recv(&(gn2[k2]), k1, MPI_DOUBLE, rightProcess, rightProcess*37, MPI_COMM_WORLD, &statusgn2);
	
	//Sending right (for sender) gna (only right part, because gna moves in x direction)
	MPI_Send(&(gna[k2-k1]),k1, MPI_DOUBLE, rightProcess, rank*41 , MPI_COMM_WORLD);
	//Recieving left (for reciever) gna (only left part, because gna moves in x direction)
	MPI_Recv(gna, k1, MPI_DOUBLE, leftProcess, leftProcess*41, MPI_COMM_WORLD, &statusgna);
	
	//Sending left (for sender) gnb (only left part, because gnb moves in -x direction)
	MPI_Send(&(gnb[k1]),k1, MPI_DOUBLE, leftProcess, rank*43 , MPI_COMM_WORLD);
	//Recieving right (for reciever) gnb (only right part, because gnb moves in -x direction)
	MPI_Recv(&(gn2[k2]), k1, MPI_DOUBLE, rightProcess, rightProcess*43, MPI_COMM_WORLD, &statusgnb);

	//Sending right (for sender) gnc (only right part, because gnc moves in x direction)
	MPI_Send(&(gnc[k2-k1]),k1, MPI_DOUBLE, rightProcess, rank*47 , MPI_COMM_WORLD);
	//Recieving left (for reciever) gn1 (only left part, because gnc moves in x direction)
	MPI_Recv(gnc, k1, MPI_DOUBLE, leftProcess, leftProcess*47, MPI_COMM_WORLD, &statusgnc);

	//Sending left (for sender) gnd (only left part, because gnd moves in -x direction)
	MPI_Send(&(gnd[k1]),k1, MPI_DOUBLE, leftProcess, rank*53 , MPI_COMM_WORLD);
	//Recieving right (for reciever) gn2 (only right part, because gnd moves in -x direction)
	MPI_Recv(&(gnd[k2]), k1, MPI_DOUBLE, rightProcess, rightProcess*53, MPI_COMM_WORLD, &statusgnd);

	//Sending right (for sender) gni (only right part, because gni moves in x direction)
	MPI_Send(&(gni[k2-k1]),k1, MPI_DOUBLE, rightProcess, rank*59 , MPI_COMM_WORLD);
	//Recieving left (for reciever) gna (only left part, because gni moves in x direction)
	MPI_Recv(gni, k1, MPI_DOUBLE, leftProcess, leftProcess*59, MPI_COMM_WORLD, &statusgni);
	
	//Sending left (for sender) gnj (only left part, because gnj moves in -x direction)
	MPI_Send(&(gnj[k1]),k1, MPI_DOUBLE, leftProcess, rank*61 , MPI_COMM_WORLD);
	//Recieving right (for reciever) gnb (only right part, because gnj moves in -x direction)
	MPI_Recv(&(gnj[k2]), k1, MPI_DOUBLE, rightProcess, rightProcess*61, MPI_COMM_WORLD, &statusgnj);

	//Sending right (for sender) gnk (only right part, because gnk moves in x direction)
	MPI_Send(&(gnk[k2-k1]),k1, MPI_DOUBLE, rightProcess, rank*67 , MPI_COMM_WORLD);
	//Recieving left (for reciever) gnk (only left part, because gnk moves in x direction)
	MPI_Recv(gnk, k1, MPI_DOUBLE, leftProcess, leftProcess*67, MPI_COMM_WORLD, &statusgnk);
	
	//Sending left (for sender) gnl (only left part, because gnl moves in -x direction)
	MPI_Send(&(gnl[k1]),k1, MPI_DOUBLE, leftProcess, rank*71 , MPI_COMM_WORLD);
	//Recieving right (for reciever) gnl (only right part, because gnl moves in -x direction)
	MPI_Recv(&(gnl[k2]), k1, MPI_DOUBLE, rightProcess, rightProcess*71, MPI_COMM_WORLD, &statusgnl);

	}
	else
	{
	
	
	//Recieving left (for reciever) fn1 (only left part, because fn1 moves in x direction)
	MPI_Recv(fn1, k1, MPI_DOUBLE, leftProcess, leftProcess*2, MPI_COMM_WORLD, &statusfn1);
	//Sending right (for sender) fn1 (only right part, because fn1 moves in x direction)
	MPI_Send(&(fn1[k2-k1]),k1, MPI_DOUBLE, rightProcess, rank*2 , MPI_COMM_WORLD);

	
	//Recieving right (for reciever) fn2 (only right part, because fn2 moves in -x direction)
	MPI_Recv(&(fn2[k2]), k1, MPI_DOUBLE, rightProcess, rightProcess*3, MPI_COMM_WORLD, &statusfn2);
	//Sending left (for sender) fn2 (only left part, because fn2 moves in -x direction)
	MPI_Send(&(fn2[k1]),k1, MPI_DOUBLE, leftProcess, rank*3 , MPI_COMM_WORLD);
	
	//Recieving left (for reciever) fna (only left part, because fna moves in x direction)
	MPI_Recv(fna, k1, MPI_DOUBLE, leftProcess, leftProcess*5, MPI_COMM_WORLD, &statusfna);
	//Sending right (for sender) fna (only right part, because fna moves in x direction)
	MPI_Send(&(fna[k2-k1]),k1, MPI_DOUBLE, rightProcess, rank*5 , MPI_COMM_WORLD);

	//Recieving right (for reciever) fnb (only right part, because fnb moves in -x direction)
	MPI_Recv(&(fn2[k2]), k1, MPI_DOUBLE, rightProcess, rightProcess*7, MPI_COMM_WORLD, &statusfnb);
	//Sending left (for sender) fnb (only left part, because fnb moves in -x direction)
	MPI_Send(&(fnb[k1]),k1, MPI_DOUBLE, leftProcess, rank*7 , MPI_COMM_WORLD);

	//Recieving left (for reciever) fn1 (only left part, because fnc moves in x direction)
	MPI_Recv(fnc, k1, MPI_DOUBLE, leftProcess, leftProcess*11, MPI_COMM_WORLD, &statusfnc);
	//Sending right (for sender) fnc (only right part, because fnc moves in x direction)
	MPI_Send(&(fnc[k2-k1]),k1, MPI_DOUBLE, rightProcess, rank*11 , MPI_COMM_WORLD);

	//Recieving right (for reciever) fn2 (only right part, because fnd moves in -x direction)
	MPI_Recv(&(fnd[k2]), k1, MPI_DOUBLE, rightProcess, rightProcess*13, MPI_COMM_WORLD, &statusfnd);
	//Sending left (for sender) fnd (only left part, because fnd moves in -x direction)
	MPI_Send(&(fnd[k1]),k1, MPI_DOUBLE, leftProcess, rank*13 , MPI_COMM_WORLD);

	//Recieving left (for reciever) fna (only left part, because fni moves in x direction)
	MPI_Recv(fni, k1, MPI_DOUBLE, leftProcess, leftProcess*17, MPI_COMM_WORLD, &statusfni);
	//Sending right (for sender) fni (only right part, because fni moves in x direction)
	MPI_Send(&(fni[k2-k1]),k1, MPI_DOUBLE, rightProcess, rank*17 , MPI_COMM_WORLD);
		
	//Recieving right (for reciever) fnb (only right part, because fnj moves in -x direction)
	MPI_Recv(&(fnj[k2]), k1, MPI_DOUBLE, rightProcess, rightProcess*19, MPI_COMM_WORLD, &statusfnj);
	//Sending left (for sender) fnj (only left part, because fnj moves in -x direction)
	MPI_Send(&(fnj[k1]),k1, MPI_DOUBLE, leftProcess, rank*19 , MPI_COMM_WORLD);

	//Recieving left (for reciever) fnk (only left part, because fnk moves in x direction)
	MPI_Recv(fnk, k1, MPI_DOUBLE, leftProcess, leftProcess*23, MPI_COMM_WORLD, &statusfnk);
	//Sending right (for sender) fnk (only right part, because fnk moves in x direction)
	MPI_Send(&(fnk[k2-k1]),k1, MPI_DOUBLE, rightProcess, rank*23 , MPI_COMM_WORLD);

	//Recieving right (for reciever) fnl (only right part, because fnl moves in -x direction)
	MPI_Recv(&(fnl[k2]), k1, MPI_DOUBLE, rightProcess, rightProcess*29, MPI_COMM_WORLD, &statusfnl);
	//Sending left (for sender) fnl (only left part, because fnl moves in -x direction)
	MPI_Send(&(fnl[k1]),k1, MPI_DOUBLE, leftProcess, rank*29 , MPI_COMM_WORLD);

	//Recieving left (for reciever) gn1 (only left part, because gn1 moves in x direction)
	MPI_Recv(gn1, k1, MPI_DOUBLE, leftProcess, leftProcess*31, MPI_COMM_WORLD, &statusgn1);
	//Sending right (for sender) gn1 (only right part, because gn1 moves in x direction)
	MPI_Send(&(gn1[k2-k1]),k1, MPI_DOUBLE, rightProcess, rank*31 , MPI_COMM_WORLD);

	//Recieving right (for reciever) gn2 (only right part, because gn2 moves in -x direction)
	MPI_Recv(&(gn2[k2]), k1, MPI_DOUBLE, rightProcess, rightProcess*37, MPI_COMM_WORLD, &statusgn2);
	//Sending left (for sender) gn2 (only left part, because gn2 moves in -x direction)
	MPI_Send(&(gn2[k1]),k1, MPI_DOUBLE, leftProcess, rank*37 , MPI_COMM_WORLD);
	
	//Recieving left (for reciever) gna (only left part, because gna moves in x direction)
	MPI_Recv(gna, k1, MPI_DOUBLE, leftProcess, leftProcess*41, MPI_COMM_WORLD, &statusgna);
	//Sending right (for sender) gna (only right part, because gna moves in x direction)
	MPI_Send(&(gna[k2-k1]),k1, MPI_DOUBLE, rightProcess, rank*41 , MPI_COMM_WORLD);
	
	//Recieving right (for reciever) gnb (only right part, because gnb moves in -x direction)
	MPI_Recv(&(gn2[k2]), k1, MPI_DOUBLE, rightProcess, rightProcess*43, MPI_COMM_WORLD, &statusgnb);
	//Sending left (for sender) gnb (only left part, because gnb moves in -x direction)
	MPI_Send(&(gnb[k1]),k1, MPI_DOUBLE, leftProcess, rank*43 , MPI_COMM_WORLD);

	//Recieving left (for reciever) gn1 (only left part, because gnc moves in x direction)
	MPI_Recv(gnc, k1, MPI_DOUBLE, leftProcess, leftProcess*47, MPI_COMM_WORLD, &statusgnc);
	//Sending right (for sender) gnc (only right part, because gnc moves in x direction)
	MPI_Send(&(gnc[k2-k1]),k1, MPI_DOUBLE, rightProcess, rank*47 , MPI_COMM_WORLD);

	//Recieving right (for reciever) gn2 (only right part, because gnd moves in -x direction)
	MPI_Recv(&(gnd[k2]), k1, MPI_DOUBLE, rightProcess, rightProcess*53, MPI_COMM_WORLD, &statusgnd);
	//Sending left (for sender) gnd (only left part, because gnd moves in -x direction)
	MPI_Send(&(gnd[k1]),k1, MPI_DOUBLE, leftProcess, rank*53 , MPI_COMM_WORLD);

	//Recieving left (for reciever) gna (only left part, because gni moves in x direction)
	MPI_Recv(gni, k1, MPI_DOUBLE, leftProcess, leftProcess*59, MPI_COMM_WORLD, &statusgni);
	//Sending right (for sender) gni (only right part, because gni moves in x direction)
	MPI_Send(&(gni[k2-k1]),k1, MPI_DOUBLE, rightProcess, rank*59 , MPI_COMM_WORLD);
	
	//Recieving right (for reciever) gnb (only right part, because gnj moves in -x direction)
	MPI_Recv(&(gnj[k2]), k1, MPI_DOUBLE, rightProcess, rightProcess*61, MPI_COMM_WORLD, &statusgnj);
	//Sending left (for sender) gnj (only left part, because gnj moves in -x direction)
	MPI_Send(&(gnj[k1]),k1, MPI_DOUBLE, leftProcess, rank*61 , MPI_COMM_WORLD);

	//Recieving left (for reciever) gnk (only left part, because gnk moves in x direction)
	MPI_Recv(gnk, k1, MPI_DOUBLE, leftProcess, leftProcess*67, MPI_COMM_WORLD, &statusgnk);
	//Sending right (for sender) gnk (only right part, because gnk moves in x direction)
	MPI_Send(&(gnk[k2-k1]),k1, MPI_DOUBLE, rightProcess, rank*67 , MPI_COMM_WORLD);
	
	//Recieving right (for reciever) gnl (only right part, because gnl moves in -x direction)
	MPI_Recv(&(gnl[k2]), k1, MPI_DOUBLE, rightProcess, rightProcess*71, MPI_COMM_WORLD, &statusgnl);
	//Sending left (for sender) gnl (only left part, because gnl moves in -x direction)
	MPI_Send(&(gnl[k1]),k1, MPI_DOUBLE, leftProcess, rank*71 , MPI_COMM_WORLD);

	}
	//if(t%infoStep==0)
		cout << "Process " << rank << ": dendities exchanged." << endl;

}

