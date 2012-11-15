#include "wet.h"

void wet::exchangeDensities(void)
{

	
	//if(t%infoStep==0)
	//	cout << "Process " << rank << ": exchanging densities...." << endl;
	
	MPI_Status statusfn1, statusfn2, statusfna, statusfnb, statusfnc, statusfnd, statusfni, statusfnj, statusfnk, statusfnl, statusgn1, statusgn2, statusgna, statusgnb, statusgnc, statusgnd, statusgni, statusgnj, statusgnk, statusgnl; 

	MPI_Request requestOutfn1, requestOutfn2, requestInfn1, requestInfn2, requestOutfna, requestOutfnb, requestInfna, requestInfnb, requestOutfnc, requestOutfnd, requestInfnc, requestInfnd, requestOutfni, requestOutfnj, requestInfni, requestInfnj, requestOutfnk, requestOutfnl, requestInfnk, requestInfnl,requestOutgn1, requestOutgn2, requestIngn1, requestIngn2, requestOutgna, requestOutgnb, requestIngna, requestIngnb, requestOutgnc, requestOutgnd, requestIngnc, requestIngnd, requestOutgni, requestOutgnj, requestIngni, requestIngnj, requestOutgnk, requestOutgnl, requestIngnk, requestIngnl;

	
	//SENDING FF DENSITIES

	//Sending right (for sender) fn1 (only right part, because fn1 moves in x direction)
	MPI_Isend(&(fn1[k2-k1]),k1, MPI_DOUBLE, rightProcess, rank*100 , MPI_COMM_WORLD, &requestOutfn1);
	
	//Sending left (for sender) fn2 (only left part, because fn2 moves in -x direction)
	MPI_Isend(&(fn2[k1]),k1, MPI_DOUBLE, leftProcess, rank*100+1 , MPI_COMM_WORLD, &requestOutfn2);

	//Sending right (for sender) fna (only right part, because fna moves in x direction)
	MPI_Isend(&(fna[k2-k1]),k1, MPI_DOUBLE, rightProcess, rank*100+2 , MPI_COMM_WORLD, &requestOutfna);
	
	//Sending left (for sender) fnb (only left part, because fnb moves in -x direction)
	MPI_Isend(&(fnb[k1]),k1, MPI_DOUBLE, leftProcess, rank*100+3 , MPI_COMM_WORLD, &requestOutfnb);

	//Sending right (for sender) fnc (only right part, because fnc moves in x direction)
	MPI_Isend(&(fnc[k2-k1]),k1, MPI_DOUBLE, rightProcess, rank*100+4 , MPI_COMM_WORLD, &requestOutfnc);
	
	//Sending left (for sender) fnd (only left part, because fnd moves in -x direction)
	MPI_Isend(&(fnd[k1]),k1, MPI_DOUBLE, leftProcess, rank*100+5 , MPI_COMM_WORLD, &requestOutfnd);

	//Sending right (for sender) fni (only right part, because fni moves in x direction)
	MPI_Isend(&(fni[k2-k1]),k1, MPI_DOUBLE, rightProcess, rank*100+6 , MPI_COMM_WORLD, &requestOutfni);
	
	//Sending left (for sender) fnj (only left part, because fnj moves in -x direction)
	MPI_Isend(&(fnj[k1]),k1, MPI_DOUBLE, leftProcess, rank*100+7 , MPI_COMM_WORLD, &requestOutfnj);

	//Sending right (for sender) fnk (only right part, because fnk moves in x direction)
	MPI_Isend(&(fnk[k2-k1]),k1, MPI_DOUBLE, rightProcess, rank*100+8 , MPI_COMM_WORLD, &requestOutfnk);
	
	//Sending left (for sender) fnl (only left part, because fnl moves in -x direction)
	MPI_Isend(&(fnl[k1]),k1, MPI_DOUBLE, leftProcess, rank*100+9 , MPI_COMM_WORLD, &requestOutfnl);

		
	//RECIEVING FF DENSITIES

	//Recieving left (for reciever) fn1 (only left part, because fn1 moves in x direction)
	MPI_Irecv(fn1, k1, MPI_DOUBLE, leftProcess, leftProcess*100, MPI_COMM_WORLD, &requestInfn1);

	//Recieving right (for reciever) fn2 (only right part, because fn2 moves in -x direction)
	MPI_Irecv(&(fn2[k2]), k1, MPI_DOUBLE, rightProcess, rightProcess*100+1, MPI_COMM_WORLD, &requestInfn2);
	
	//Recieving left (for reciever) fna (only left part, because fna moves in x direction)
	MPI_Irecv(fna, k1, MPI_DOUBLE, leftProcess, leftProcess*100+2, MPI_COMM_WORLD, &requestInfna);
	
	//Recieving right (for reciever) fnb (only right part, because fnb moves in -x direction)
	MPI_Irecv(&(fnb[k2]), k1, MPI_DOUBLE, rightProcess, rightProcess*100+3, MPI_COMM_WORLD, &requestInfnb);

	//Recieving left (for reciever) fn1 (only left part, because fnc moves in x direction)
	MPI_Irecv(fnc, k1, MPI_DOUBLE, leftProcess, leftProcess*100+4, MPI_COMM_WORLD, &requestInfnc);

	//Recieving right (for reciever) fn2 (only right part, because fnd moves in -x direction)
	MPI_Irecv(&(fnd[k2]), k1, MPI_DOUBLE, rightProcess, rightProcess*100+5, MPI_COMM_WORLD, &requestInfnd);
	
	//Recieving left (for reciever) fna (only left part, because fni moves in x direction)
	MPI_Irecv(fni, k1, MPI_DOUBLE, leftProcess, leftProcess*100+6, MPI_COMM_WORLD, &requestInfni);
	
	//Recieving right (for reciever) fnb (only right part, because fnj moves in -x direction)
	MPI_Irecv(&(fnj[k2]), k1, MPI_DOUBLE, rightProcess, rightProcess*100+7, MPI_COMM_WORLD, &requestInfnj);

	//Recieving left (for reciever) fnk (only left part, because fnk moves in x direction)
	MPI_Irecv(fnk, k1, MPI_DOUBLE, leftProcess, leftProcess*100+8, MPI_COMM_WORLD, &requestInfnk);
	
	//Recieving right (for reciever) fnl (only right part, because fnl moves in -x direction)
	MPI_Irecv(&(fnl[k2]), k1, MPI_DOUBLE, rightProcess, rightProcess*100+9, MPI_COMM_WORLD, &requestInfnl);

	
	//SENDING GG DENSITIES

	//Sending right (for sender) gn1 (only right part, because gn1 moves in x direction)
	MPI_Isend(&(gn1[k2-k1]),k1, MPI_DOUBLE, rightProcess, rank*100+10 , MPI_COMM_WORLD, &requestOutgn1);
	
	//Sending left (for sender) gn2 (only left part, because gn2 moves in -x direction)
	MPI_Isend(&(gn2[k1]),k1, MPI_DOUBLE, leftProcess, rank*100+11 , MPI_COMM_WORLD, &requestOutgn2);

	//Sending right (for sender) gna (only right part, because gna moves in x direction)
	MPI_Isend(&(gna[k2-k1]),k1, MPI_DOUBLE, rightProcess, rank*100+12 , MPI_COMM_WORLD, &requestOutgna);
	
	//Sending left (for sender) gnb (only left part, because gnb moves in -x direction)
	MPI_Isend(&(gnb[k1]),k1, MPI_DOUBLE, leftProcess, rank*100+13 , MPI_COMM_WORLD, &requestOutgnb);

	//Sending right (for sender) gnc (only right part, because gnc moves in x direction)
	MPI_Isend(&(gnc[k2-k1]),k1, MPI_DOUBLE, rightProcess, rank*100+14 , MPI_COMM_WORLD, &requestOutgnc);
	
	//Sending left (for sender) gnd (only left part, because gnd moves in -x direction)
	MPI_Isend(&(gnd[k1]),k1, MPI_DOUBLE, leftProcess, rank*100+15, MPI_COMM_WORLD, &requestOutgnd);

	//Sending right (for sender) gni (only right part, because gni moves in x direction)
	MPI_Isend(&(gni[k2-k1]),k1, MPI_DOUBLE, rightProcess, rank*100+16 , MPI_COMM_WORLD, &requestOutgni);
	
	//Sending left (for sender) gnj (only left part, because gnj moves in -x direction)
	MPI_Isend(&(gnj[k1]),k1, MPI_DOUBLE, leftProcess, rank*100+17 , MPI_COMM_WORLD, &requestOutgnj);

	//Sending right (for sender) gnk (only right part, because gnk moves in x direction)
	MPI_Isend(&(gnk[k2-k1]),k1, MPI_DOUBLE, rightProcess, rank*100+18 , MPI_COMM_WORLD, &requestOutgnk);
	
	//Sending left (for sender) gnl (only left part, because gnl moves in -x direction)
	MPI_Isend(&(gnl[k1]),k1, MPI_DOUBLE, leftProcess, rank*100+19 , MPI_COMM_WORLD, &requestOutgnl);

		
	//RECIEVING GG DENSITIES

	//Recieving left (for reciever) gn1 (only left part, because gn1 moves in x direction)
	MPI_Irecv(gn1, k1, MPI_DOUBLE, leftProcess, leftProcess*100+10, MPI_COMM_WORLD, &requestIngn1);

	//Recieving right (for reciever) gn2 (only right part, because gn2 moves in -x direction)
	MPI_Irecv(&(gn2[k2]), k1, MPI_DOUBLE, rightProcess, rightProcess*100+11, MPI_COMM_WORLD, &requestIngn2);
	
	//Recieving left (for reciever) gna (only left part, because gna moves in x direction)
	MPI_Irecv(gna, k1, MPI_DOUBLE, leftProcess, leftProcess*100+12, MPI_COMM_WORLD, &requestIngna);
	
	//Recieving right (for reciever) gnb (only right part, because gnb moves in -x direction)
	MPI_Irecv(&(gnb[k2]), k1, MPI_DOUBLE, rightProcess, rightProcess*100+13, MPI_COMM_WORLD, &requestIngnb);

	//Recieving left (for reciever) gn1 (only left part, because gnc moves in x direction)
	MPI_Irecv(gnc, k1, MPI_DOUBLE, leftProcess, leftProcess*100+14, MPI_COMM_WORLD, &requestIngnc);

	//Recieving right (for reciever) gn2 (only right part, because gnd moves in -x direction)
	MPI_Irecv(&(gnd[k2]), k1, MPI_DOUBLE, rightProcess, rightProcess*100+15, MPI_COMM_WORLD, &requestIngnd);
	
	//Recieving left (for reciever) gna (only left part, because gni moves in x direction)
	MPI_Irecv(gni, k1, MPI_DOUBLE, leftProcess, leftProcess*100+16, MPI_COMM_WORLD, &requestIngni);
	
	//Recieving right (for reciever) gnb (only right part, because gnj moves in -x direction)
	MPI_Irecv(&(gnj[k2]), k1, MPI_DOUBLE, rightProcess, rightProcess*100+17, MPI_COMM_WORLD, &requestIngnj);

	//Recieving left (for reciever) gnk (only left part, because gnk moves in x direction)
	MPI_Irecv(gnk, k1, MPI_DOUBLE, leftProcess, leftProcess*100+18, MPI_COMM_WORLD, &requestIngnk);
	
	//Recieving right (for reciever) gnl (only right part, because gnl moves in -x direction)
	MPI_Irecv(&(gnl[k2]), k1, MPI_DOUBLE, rightProcess, rightProcess*100+19, MPI_COMM_WORLD, &requestIngnl);

	
	
	//Waiting until fn(i) densities are saved in the recieving bufner
	MPI_Wait(&requestInfn1, &statusfn1);
	MPI_Wait(&requestInfn2, &statusfn2);
	MPI_Wait(&requestInfna, &statusfna);
	MPI_Wait(&requestInfnb, &statusfnb);
	MPI_Wait(&requestInfnc, &statusfnc);
	MPI_Wait(&requestInfnd, &statusfnd);
	MPI_Wait(&requestInfni, &statusfni);
	MPI_Wait(&requestInfnj, &statusfnj);
	MPI_Wait(&requestInfnk, &statusfnk);
	MPI_Wait(&requestInfnl, &statusfnl);

	//Waiting until gn(i) densities are saved in the recieving bufner
	MPI_Wait(&requestIngn1, &statusgn1);
	MPI_Wait(&requestIngn2, &statusgn2);
	MPI_Wait(&requestIngna, &statusgna);
	MPI_Wait(&requestIngnb, &statusgnb);
	MPI_Wait(&requestIngnc, &statusgnc);
	MPI_Wait(&requestIngnd, &statusgnd);
	MPI_Wait(&requestIngni, &statusgni);
	MPI_Wait(&requestIngnj, &statusgnj);
	MPI_Wait(&requestIngnk, &statusgnk);
	MPI_Wait(&requestIngnl, &statusgnl);

	//Waiting until fn(i) sending bufner is relased
	MPI_Wait(&requestOutfn1, &statusfn1);
	MPI_Wait(&requestOutfn2, &statusfn2);
	MPI_Wait(&requestOutfna, &statusfna);
	MPI_Wait(&requestOutfnb, &statusfnb);
	MPI_Wait(&requestOutfnc, &statusfnc);
	MPI_Wait(&requestOutfnd, &statusfnd);
	MPI_Wait(&requestOutfni, &statusfni);
	MPI_Wait(&requestOutfnj, &statusfnj);
	MPI_Wait(&requestOutfnk, &statusfnk);
	MPI_Wait(&requestOutfnl, &statusfnl);

	//Waiting until gn(i) sending bufner is relased
	MPI_Wait(&requestOutgn1, &statusgn1);
	MPI_Wait(&requestOutgn2, &statusgn2);
	MPI_Wait(&requestOutgna, &statusgna);
	MPI_Wait(&requestOutgnb, &statusgnb);
	MPI_Wait(&requestOutgnc, &statusgnc);
	MPI_Wait(&requestOutgnd, &statusgnd);
	MPI_Wait(&requestOutgni, &statusgni);
	MPI_Wait(&requestOutgnj, &statusgnj);
	MPI_Wait(&requestOutgnk, &statusgnk);
	MPI_Wait(&requestOutgnl, &statusgnl);

	MPI_Barrier(MPI_COMM_WORLD);

	//if(t%infoStep==0)
	//	cout << "Process " << rank << ": dendities exchanged." << endl;
}

