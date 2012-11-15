//Costructor: initialize variables defined in the class wet
//Marco, 29 June 2011, Lisa, 28 July 2011

#include "wet.h"


wet::wet(void)
{
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);		//rank = id number of process (see mpi.h library)
	MPI_Comm_size(MPI_COMM_WORLD, &size);		//size = total number of processes (see mpi.h library)

	cout << "Process " << rank << ": Constructor: initialisation..." << endl;

	t=0;				//index 				(SISTEMA)
	afterequilflag=0;
	
	//put into initialise, Lisa 31/07/2011, for equilibration/duplication
	/*
	if(rank==0)
		leftProcess=size-1;
	else
		leftProcess=rank-1;

	if(rank==size-1)
		rightProcess=0;
	else
		rightProcess=rank+1;
	*/
	//cout << "Constructor: Process " << rank << " : rightProcess: " << rightProcess << " leftProcess: " << leftProcess << endl;
	ReadInput();
    
   
	
    
    if(rank==ROOT)
	{
		
        char command[50];
        snprintf(command, 50,"mkdir uxs%.3frat%.3foff%.0f", initUX, dropletR/Dx, dropletCenterY-Px-Dx/2);
        
        cout <<"String command is equal to " << command << endl;
        
        
        //system("rm -r densitydata");
		system(command);		//create a directory where store data
		//system("rm *.dat");
		

		//ReadInput();			//process ROOT ead data from wet.dat (lattice structure, number of steps, ...), from rectpost.at (solid surface structure), and from spherical.dat (droplet characteristic)
	}

	
	//commReadData();		//comunicate data read by process ROOT
	
	initialise();			//Initialise variables, surface and gas configuration;
	
	
	//Calculate the surface free energy constants 
	//Here phi11 = - phi1 / kappa (in literarure (Kusumaatmaja, Yeomans, Lattice Boltzmann Simulations of Wetting and Drop Dynamics)) 
	//And  phi12 = - phi2 / kappa (in literature (Kusumaatmaja, Yeomans, Lattice Boltzmann Simulations of Wetting and Drop Dynamics))
	
	

	
/*
	energy=0.0;			//energy, used in LBAlgorithm()						(SISTEMA)
	bulkE=0.0;			//bulk energy, used in LBAlgorithm()					(SISTEMA)
	interfaceE=0.0;			//interface energy, used in LBAlgorithm()				(SISTEMA)
	surfaceE=0.0;			//surface energy, used in LBAlgorithm()					(SISTEMA)
	//these two are just for checking whether my (Lisa) freeE calculation is correct			(SISTEMA)
	intE_5h=0.0; intE_3h=0.0; intE_5H=0.0; intE_10H=0.0;    
	//used in LBAlgorithm()											(SISTEMA)
	gradmax=0.0; gradmax2=0.0; gradmaxx=0.0; gradmax2x=0.0; gradmaxy=0.0; gradmax2y=0.0; gradmaxz=0.0; gradmax2z=0.0; 
	
	
	RXcm=0.0, RXcmOld=0.0;		//Center of mass of fluid f and the same "old" variable along x, in ARolling()
	RZcm=0.0;			//Center of mass of fluid f along z, in ARolling()
	RNtot=0.0;			//TOTAL mass, sum of n for all N lattice points, calculated in ARolling()
	RVXcm=0.0; RVYcm=0.0; RVZcm=0.0;//Velocity of center of mass of fluid f, in ARolling()
	kin=0.0; 			//TOTAL kinetic energy of fluid f, calculated in ARolling()
	kin_CM=0.0;			//kinetic energy of center of mass of fluid f, calculated in ARolling()	
	ratio_kin=0.0;			//ratio between kinetic energy of center of mass and total kinetic energy (kin_CM/kin), in ARolling()
	d_xk=0.0; d_yk=0.0; d_zk=0.0;	//lattice points, 3-D index version of index k (like xk, yk, zk, but double), in ARolling()
*/
	cout << "Process " << rank << ": Constructor: data initialized" << endl;
}

