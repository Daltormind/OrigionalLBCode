#include "wet.h"

void wet::exchangeDensities_ffgg(void)
{

	
	//if(t%infoStep==0)
	//	cout << "Process " << rank << ": exchanging densities...." << endl;
	
	MPI_Status statusff1, statusff2, statusffa, statusffb, statusffc, statusffd, statusffi, statusffj, statusffk, statusffl, statusgg1, statusgg2, statusgga, statusggb, statusggc, statusggd, statusggi, statusggj, statusggk, statusggl; 

	MPI_Request requestOutff1, requestOutff2, requestInff1, requestInff2, requestOutffa, requestOutffb, requestInffa, requestInffb, requestOutffc, requestOutffd, requestInffc, requestInffd, requestOutffi, requestOutffj, requestInffi, requestInffj, requestOutffk, requestOutffl, requestInffk, requestInffl,requestOutgg1, requestOutgg2, requestIngg1, requestIngg2, requestOutgga, requestOutggb, requestIngga, requestInggb, requestOutggc, requestOutggd, requestInggc, requestInggd, requestOutggi, requestOutggj, requestInggi, requestInggj, requestOutggk, requestOutggl, requestInggk, requestInggl;

	
	//SENDING FF DENSITIES

	//Sending right (for sender) ff1 (only right part, because ff1 moves in x direction)
	MPI_Isend(&(ff1[k2-k1]),k1, MPI_DOUBLE, rightProcess, rank*100 , MPI_COMM_WORLD, &requestOutff1);
	
	//Sending left (for sender) ff2 (only left part, because ff2 moves in -x direction)
	MPI_Isend(&(ff2[k1]),k1, MPI_DOUBLE, leftProcess, rank*100+1 , MPI_COMM_WORLD, &requestOutff2);

	//Sending right (for sender) ffa (only right part, because ffa moves in x direction)
	MPI_Isend(&(ffa[k2-k1]),k1, MPI_DOUBLE, rightProcess, rank*100+2 , MPI_COMM_WORLD, &requestOutffa);
	
	//Sending left (for sender) ffb (only left part, because ffb moves in -x direction)
	MPI_Isend(&(ffb[k1]),k1, MPI_DOUBLE, leftProcess, rank*100+3 , MPI_COMM_WORLD, &requestOutffb);

	//Sending right (for sender) ffc (only right part, because ffc moves in x direction)
	MPI_Isend(&(ffc[k2-k1]),k1, MPI_DOUBLE, rightProcess, rank*100+4 , MPI_COMM_WORLD, &requestOutffc);
	
	//Sending left (for sender) ffd (only left part, because ffd moves in -x direction)
	MPI_Isend(&(ffd[k1]),k1, MPI_DOUBLE, leftProcess, rank*100+5 , MPI_COMM_WORLD, &requestOutffd);

	//Sending right (for sender) ffi (only right part, because ffi moves in x direction)
	MPI_Isend(&(ffi[k2-k1]),k1, MPI_DOUBLE, rightProcess, rank*100+6 , MPI_COMM_WORLD, &requestOutffi);
	
	//Sending left (for sender) ffj (only left part, because ffj moves in -x direction)
	MPI_Isend(&(ffj[k1]),k1, MPI_DOUBLE, leftProcess, rank*100+7 , MPI_COMM_WORLD, &requestOutffj);

	//Sending right (for sender) ffk (only right part, because ffk moves in x direction)
	MPI_Isend(&(ffk[k2-k1]),k1, MPI_DOUBLE, rightProcess, rank*100+8 , MPI_COMM_WORLD, &requestOutffk);
	
	//Sending left (for sender) ffl (only left part, because ffl moves in -x direction)
	MPI_Isend(&(ffl[k1]),k1, MPI_DOUBLE, leftProcess, rank*100+9 , MPI_COMM_WORLD, &requestOutffl);

		
	//RECIEVING FF DENSITIES

	//Recieving left (for reciever) ff1 (only left part, because ff1 moves in x direction)
	MPI_Irecv(ff1, k1, MPI_DOUBLE, leftProcess, leftProcess*100, MPI_COMM_WORLD, &requestInff1);

	//Recieving right (for reciever) ff2 (only right part, because ff2 moves in -x direction)
	MPI_Irecv(&(ff2[k2]), k1, MPI_DOUBLE, rightProcess, rightProcess*100+1, MPI_COMM_WORLD, &requestInff2);
	
	//Recieving left (for reciever) ffa (only left part, because ffa moves in x direction)
	MPI_Irecv(ffa, k1, MPI_DOUBLE, leftProcess, leftProcess*100+2, MPI_COMM_WORLD, &requestInffa);
	
	//Recieving right (for reciever) ffb (only right part, because ffb moves in -x direction)
	MPI_Irecv(&(ffb[k2]), k1, MPI_DOUBLE, rightProcess, rightProcess*100+3, MPI_COMM_WORLD, &requestInffb);

	//Recieving left (for reciever) ff1 (only left part, because ffc moves in x direction)
	MPI_Irecv(ffc, k1, MPI_DOUBLE, leftProcess, leftProcess*100+4, MPI_COMM_WORLD, &requestInffc);

	//Recieving right (for reciever) ff2 (only right part, because ffd moves in -x direction)
	MPI_Irecv(&(ffd[k2]), k1, MPI_DOUBLE, rightProcess, rightProcess*100+5, MPI_COMM_WORLD, &requestInffd);
	
	//Recieving left (for reciever) ffa (only left part, because ffi moves in x direction)
	MPI_Irecv(ffi, k1, MPI_DOUBLE, leftProcess, leftProcess*100+6, MPI_COMM_WORLD, &requestInffi);
	
	//Recieving right (for reciever) ffb (only right part, because ffj moves in -x direction)
	MPI_Irecv(&(ffj[k2]), k1, MPI_DOUBLE, rightProcess, rightProcess*100+7, MPI_COMM_WORLD, &requestInffj);

	//Recieving left (for reciever) ffk (only left part, because ffk moves in x direction)
	MPI_Irecv(ffk, k1, MPI_DOUBLE, leftProcess, leftProcess*100+8, MPI_COMM_WORLD, &requestInffk);
	
	//Recieving right (for reciever) ffl (only right part, because ffl moves in -x direction)
	MPI_Irecv(&(ffl[k2]), k1, MPI_DOUBLE, rightProcess, rightProcess*100+9, MPI_COMM_WORLD, &requestInffl);

	
	//SENDING GG DENSITIES

	//Sending right (for sender) gg1 (only right part, because gg1 moves in x direction)
	MPI_Isend(&(gg1[k2-k1]),k1, MPI_DOUBLE, rightProcess, rank*100+10 , MPI_COMM_WORLD, &requestOutgg1);
	
	//Sending left (for sender) gg2 (only left part, because gg2 moves in -x direction)
	MPI_Isend(&(gg2[k1]),k1, MPI_DOUBLE, leftProcess, rank*100+11 , MPI_COMM_WORLD, &requestOutgg2);

	//Sending right (for sender) gga (only right part, because gga moves in x direction)
	MPI_Isend(&(gga[k2-k1]),k1, MPI_DOUBLE, rightProcess, rank*100+12 , MPI_COMM_WORLD, &requestOutgga);
	
	//Sending left (for sender) ggb (only left part, because ggb moves in -x direction)
	MPI_Isend(&(ggb[k1]),k1, MPI_DOUBLE, leftProcess, rank*100+13 , MPI_COMM_WORLD, &requestOutggb);

	//Sending right (for sender) ggc (only right part, because ggc moves in x direction)
	MPI_Isend(&(ggc[k2-k1]),k1, MPI_DOUBLE, rightProcess, rank*100+14 , MPI_COMM_WORLD, &requestOutggc);
	
	//Sending left (for sender) ggd (only left part, because ggd moves in -x direction)
	MPI_Isend(&(ggd[k1]),k1, MPI_DOUBLE, leftProcess, rank*100+15, MPI_COMM_WORLD, &requestOutggd);

	//Sending right (for sender) ggi (only right part, because ggi moves in x direction)
	MPI_Isend(&(ggi[k2-k1]),k1, MPI_DOUBLE, rightProcess, rank*100+16 , MPI_COMM_WORLD, &requestOutggi);
	
	//Sending left (for sender) ggj (only left part, because ggj moves in -x direction)
	MPI_Isend(&(ggj[k1]),k1, MPI_DOUBLE, leftProcess, rank*100+17 , MPI_COMM_WORLD, &requestOutggj);

	//Sending right (for sender) ggk (only right part, because ggk moves in x direction)
	MPI_Isend(&(ggk[k2-k1]),k1, MPI_DOUBLE, rightProcess, rank*100+18 , MPI_COMM_WORLD, &requestOutggk);
	
	//Sending left (for sender) ggl (only left part, because ggl moves in -x direction)
	MPI_Isend(&(ggl[k1]),k1, MPI_DOUBLE, leftProcess, rank*100+19 , MPI_COMM_WORLD, &requestOutggl);

		
	//RECIEVING GG DENSITIES

	//Recieving left (for reciever) gg1 (only left part, because gg1 moves in x direction)
	MPI_Irecv(gg1, k1, MPI_DOUBLE, leftProcess, leftProcess*100+10, MPI_COMM_WORLD, &requestIngg1);

	//Recieving right (for reciever) gg2 (only right part, because gg2 moves in -x direction)
	MPI_Irecv(&(gg2[k2]), k1, MPI_DOUBLE, rightProcess, rightProcess*100+11, MPI_COMM_WORLD, &requestIngg2);
	
	//Recieving left (for reciever) gga (only left part, because gga moves in x direction)
	MPI_Irecv(gga, k1, MPI_DOUBLE, leftProcess, leftProcess*100+12, MPI_COMM_WORLD, &requestIngga);
	
	//Recieving right (for reciever) ggb (only right part, because ggb moves in -x direction)
	MPI_Irecv(&(ggb[k2]), k1, MPI_DOUBLE, rightProcess, rightProcess*100+13, MPI_COMM_WORLD, &requestInggb);

	//Recieving left (for reciever) gg1 (only left part, because ggc moves in x direction)
	MPI_Irecv(ggc, k1, MPI_DOUBLE, leftProcess, leftProcess*100+14, MPI_COMM_WORLD, &requestInggc);

	//Recieving right (for reciever) gg2 (only right part, because ggd moves in -x direction)
	MPI_Irecv(&(ggd[k2]), k1, MPI_DOUBLE, rightProcess, rightProcess*100+15, MPI_COMM_WORLD, &requestInggd);
	
	//Recieving left (for reciever) gga (only left part, because ggi moves in x direction)
	MPI_Irecv(ggi, k1, MPI_DOUBLE, leftProcess, leftProcess*100+16, MPI_COMM_WORLD, &requestInggi);
	
	//Recieving right (for reciever) ggb (only right part, because ggj moves in -x direction)
	MPI_Irecv(&(ggj[k2]), k1, MPI_DOUBLE, rightProcess, rightProcess*100+17, MPI_COMM_WORLD, &requestInggj);

	//Recieving left (for reciever) ggk (only left part, because ggk moves in x direction)
	MPI_Irecv(ggk, k1, MPI_DOUBLE, leftProcess, leftProcess*100+18, MPI_COMM_WORLD, &requestInggk);
	
	//Recieving right (for reciever) ggl (only right part, because ggl moves in -x direction)
	MPI_Irecv(&(ggl[k2]), k1, MPI_DOUBLE, rightProcess, rightProcess*100+19, MPI_COMM_WORLD, &requestInggl);

	
	
	//Waiting until ff(i) densities are saved in the recieving buffer
	MPI_Wait(&requestInff1, &statusff1);
	MPI_Wait(&requestInff2, &statusff2);
	MPI_Wait(&requestInffa, &statusffa);
	MPI_Wait(&requestInffb, &statusffb);
	MPI_Wait(&requestInffc, &statusffc);
	MPI_Wait(&requestInffd, &statusffd);
	MPI_Wait(&requestInffi, &statusffi);
	MPI_Wait(&requestInffj, &statusffj);
	MPI_Wait(&requestInffk, &statusffk);
	MPI_Wait(&requestInffl, &statusffl);

	//Waiting until gg(i) densities are saved in the recieving buffer
	MPI_Wait(&requestIngg1, &statusgg1);
	MPI_Wait(&requestIngg2, &statusgg2);
	MPI_Wait(&requestIngga, &statusgga);
	MPI_Wait(&requestInggb, &statusggb);
	MPI_Wait(&requestInggc, &statusggc);
	MPI_Wait(&requestInggd, &statusggd);
	MPI_Wait(&requestInggi, &statusggi);
	MPI_Wait(&requestInggj, &statusggj);
	MPI_Wait(&requestInggk, &statusggk);
	MPI_Wait(&requestInggl, &statusggl);

	//Waiting until ff(i) sending buffer is relased
	MPI_Wait(&requestOutff1, &statusff1);
	MPI_Wait(&requestOutff2, &statusff2);
	MPI_Wait(&requestOutffa, &statusffa);
	MPI_Wait(&requestOutffb, &statusffb);
	MPI_Wait(&requestOutffc, &statusffc);
	MPI_Wait(&requestOutffd, &statusffd);
	MPI_Wait(&requestOutffi, &statusffi);
	MPI_Wait(&requestOutffj, &statusffj);
	MPI_Wait(&requestOutffk, &statusffk);
	MPI_Wait(&requestOutffl, &statusffl);

	//Waiting until gg(i) sending buffer is relased
	MPI_Wait(&requestOutgg1, &statusgg1);
	MPI_Wait(&requestOutgg2, &statusgg2);
	MPI_Wait(&requestOutgga, &statusgga);
	MPI_Wait(&requestOutggb, &statusggb);
	MPI_Wait(&requestOutggc, &statusggc);
	MPI_Wait(&requestOutggd, &statusggd);
	MPI_Wait(&requestOutggi, &statusggi);
	MPI_Wait(&requestOutggj, &statusggj);
	MPI_Wait(&requestOutggk, &statusggk);
	MPI_Wait(&requestOutggl, &statusggl);

	MPI_Barrier(MPI_COMM_WORLD);

	//if(t%infoStep==0)
	//	cout << "Process " << rank << ": dendities exchanged." << endl;
}

