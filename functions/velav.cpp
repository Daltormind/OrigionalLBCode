//velav works out average velocity of fluid

#include "wet.h"

void wet::velav()
{

double sumx,sumy,sumz,num,tx,ty,tz,thresh;


sumx=0.0;
sumy=0.0;
sumz=0.0;
num=0.0;

MPI_Status statusROOT;

MPI_Gather(&(uxs[k1]),k2-k1,MPI_DOUBLE,uGlobal,k2-k1,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);

if(rank==ROOT)
{
for( k = 0 ; k < N ; k++)
{         
	if(pGlobal[k]>p_thresh_n)
	{
		sumx+=uGlobal[k];
		num+=1;
	}
}

avx=sumx/num;
}

MPI_Gather(&(uys[k1]),k2-k1,MPI_DOUBLE,uGlobal,k2-k1,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);

if(rank==ROOT)
{
for( k = 0 ; k < N ; k++)
{         
	if(pGlobal[k]>p_thresh_n)
	{
		sumy+=uGlobal[k];
		
	}
}



avy=sumy/num;
}

MPI_Gather(&(uzs[k1]),k2-k1,MPI_DOUBLE,uGlobal,k2-k1,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);

if(rank==ROOT)
{
for( k = 0 ; k < N ; k++)
{         
	if(pGlobal[k]>p_thresh_n)
	{
		sumz+=uGlobal[k];
		num+=1;
	}
}

avz=sumz/num;

tx=abs(avx-avxt);
ty=abs(avy-avyt);
tz=abs(avz-avzt);

if(t%writeStep==0)
{
	cout << "t " << t << " tx ty tz " << tx << " "<< ty << " "<< tz << endl;
}

thresh=0.000001;
if(tx<thresh and ty<thresh and tz<thresh)
{
	co+=1;
	//cout << "co= " << co << endl;
}
else
{
	avxt=avx;
	avyt=avy;
	avzt=avz;
	co=0;
}
}

if(rank==ROOT)
{
	for(int r=1;r<size;r++)
	{
		MPI_Send(&co, 1, MPI_INT, r, rank, MPI_COMM_WORLD);
	}
}
else
{
	MPI_Recv(&co, 1, MPI_INT, ROOT, ROOT, MPI_COMM_WORLD,&statusROOT); 
}



}