//main.cpp: main file to execute the simulation
//Marco, 29 June 2011

#include "wet.h"
#include "mpi.h"
#include <iostream>
//#include "String.hh"
//#include "Bool.hh"

using namespace std;


int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);
	
	wet drop;
	cout << " --- Starting LB Algorithm ---" << endl;
	
		
		//cout << "--- Process " << drop.rank <<": nbEqStep "<< drop.nbEqStep << endl;
		drop.LBAlgorithm(0,drop.nbEqStep);
	

	
	MPI_Finalize();

	return 0;
}

