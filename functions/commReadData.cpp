#include "wet.h"
#include <iostream>

using namespace std;

void wet::commReadData(void)
{

	if(MPI_Bcast(&equilfirst, 1, MPI_INT, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS)
		cout << "Process "<< rank << ": equilfirst = " << equilfirst << endl;	
	if(MPI_Bcast(&velocityprofile, 1, MPI_INT, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS)
		cout << "Process "<< rank << ": velocityprofile = " << velocityprofile << endl;
	if(MPI_Bcast(&duplicationtype, 1, MPI_INT, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS)
		cout << "Process "<< rank << ": duplicationtype = " << duplicationtype << endl;
	if(MPI_Bcast(&dimensions, 1, MPI_INT, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS)
		cout << "Process "<< rank << ": dimensions = " << dimensions << endl;
	if(MPI_Bcast(&LXc, 1, MPI_INT, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS)
		cout << "Process "<< rank << ": LXc = " << LXc << endl;
	if(MPI_Bcast(&LY, 1, MPI_INT, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS)
		cout << "Process "<< rank << ": LY = " << LY << endl;
	if(MPI_Bcast(&LZ, 1, MPI_INT, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS)
		cout << "Process "<< rank << ": LZ = " << LZ << endl;
	if(MPI_Bcast(&nbEqStep, 1, MPI_INT, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS)
		cout << "Process "<< rank << ": nbEqStep = " << nbEqStep << endl;
	if(MPI_Bcast(&equilTime, 1, MPI_INT, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS)
		cout << "Process "<< rank << ": equilTime = " << equilTime << endl;
	if(MPI_Bcast(&infoStep, 1, MPI_INT, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS)
		cout << "Process "<< rank << ": infoStep = " << infoStep << endl;
	if(MPI_Bcast(&writeStep, 1, MPI_INT, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS)
		cout << "Process "<< rank << ": writeStep = " << writeStep << endl; 
	if(MPI_Bcast(&floorTime, 1, MPI_INT, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS)
		cout << "Process "<< rank << ": floorTime = " << floorTime << endl;
	if(MPI_Bcast(&teta1, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS)
		cout << "Process "<< rank << ": teta1 = " << teta1 << endl;
	if(MPI_Bcast(&teta2, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS)
		cout << "Process "<< rank << ": teta2 = " << teta2 << endl;
	if(MPI_Bcast(&teta_CB_fix, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS)
		cout << "Process "<< rank << ": teta_CB_fix = " << teta_CB_fix << endl;
	if(MPI_Bcast(&tauliquid, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS)
		cout << "Process "<< rank << ": tauliquid = " << tauliquid << endl;
	if(MPI_Bcast(&taugas, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS)
		cout << "Process "<< rank << ": taugas = " << taugas << endl;
	if(MPI_Bcast(&tau1, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS)
		cout << "Process "<< rank << ": tau1 = " << tau1 << endl;
	if(MPI_Bcast(&kappa_p, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS)
		cout << "Process "<< rank << ": kappa_p = " << kappa_p << endl;
	if(MPI_Bcast(&A, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS)
		cout << "Process "<< rank << ": A = " << A << endl;
	if(MPI_Bcast(&B, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS)
		cout << "Process "<< rank << ": B = " << B << endl;
	if(MPI_Bcast(&gama, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS)
		cout << "Process "<< rank << ": gama = " << gama << endl;
	if(MPI_Bcast(G, 3, MPI_DOUBLE, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS)
		cout << "Process "<< rank << ": F_x = " << G[0] << " F_y = " << G[1] << " F_z = " << G[2] << endl;
	if(MPI_Bcast(&p_thresh, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS)
		cout << "Process "<< rank << ": p_thresh = " << p_thresh << endl;
		
	
	if(MPI_Bcast(&geometry, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS)
		cout << "Process "<< rank << ": geometry = " << geometry << endl;
	if(MPI_Bcast(&Dx, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS)
		cout << "Process "<< rank << ": Dx = " << Dx << endl;
	if(MPI_Bcast(&Dy, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS)
		cout << "Process "<< rank << ": Dy = " << Dy << endl;
	if(MPI_Bcast(&Dh, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS)
		cout << "Process "<< rank << ": Dh = " << Dh << endl;
	if(MPI_Bcast(&PeriodX, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS)
		cout << "Process "<< rank << ": PeriodX = " << PeriodX << endl;
	if(MPI_Bcast(&PeriodY, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS)
		cout << "Process "<< rank << ": PeriodY = " << PeriodY << endl;
	if(MPI_Bcast(&PostSymm, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS)
		cout << "Process "<< rank << ": PostSymm = " << tau1 << endl;



	if(MPI_Bcast(&dropletR, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS)
		cout << "Process "<< rank << ": dropletR = " << dropletR << endl;
	if(MPI_Bcast(&dropletCenterX, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS)
		cout << "Process "<< rank << ": dropletCenterX = " << dropletCenterX << endl;
	if(MPI_Bcast(&dropletCenterY, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS)
		cout << "Process "<< rank << ": dropletCenterY = " << dropletCenterY << endl;
	if(MPI_Bcast(&dropletCenterZ, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS)
		cout << "Process "<< rank << ": dropletCenterZ = " << dropletCenterZ << endl;
	if(MPI_Bcast(&initUX, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS)
		cout << "Process "<< rank << ": initUX = " << initUX << endl;
	if(MPI_Bcast(&initUY, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS)
		cout << "Process "<< rank << ": initUY = " << initUY << endl;
	if(MPI_Bcast(&initUZ, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD) == MPI_SUCCESS)
		cout << "Process "<< rank << ": initUZ = " << initUZ << endl;

}
