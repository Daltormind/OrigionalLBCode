//ReadInput.cpp: read input from wet.dat, rectpost.dat and spherical.dat
//Marco, 28 June 2011, Lisa, 28 July 2011

#include "wet.h"


int wet::ReadInput(void)
{	
	
	String endOfLine; 
	
	ifstream inputFile0("settings.par");
	if (!inputFile0) {
		cout << "Can't open the input file settings.par " << endl;
	}
	inputFile0 >> equilfirst >> endOfLine;		if (rank==ROOT) cout << "equilfirst? " << equilfirst << endl;
	inputFile0 >> velocityprofile >> endOfLine;	if (rank==ROOT) cout << "velocityprofile? " << velocityprofile << endl;
	inputFile0 >> duplicationtype  >> endOfLine;	if (rank==ROOT) cout << "duplicationtype " << duplicationtype  << endl;
	inputFile0.close();	 

	ifstream inputFile("wet.par");
	if (!inputFile) {
		cout << "Can't open the input file wet.par" << endl;
		return 1;
	}
	
	inputFile >> dimensions >> endOfLine;		if (rank==ROOT) cout << "dimensions " << dimensions << endl;  
	if (dimensions != 2 && dimensions !=3) {
		cout << "Dimensions not well defined in settings.par (2 or 3?)" << endl;
		return 1;
	}
	inputFile >> LXc >> endOfLine;				if (rank==ROOT) cout << "LX " << LXc << endl;   //Lisa: LXc is the constant value that is not to be changed. LX is changed for equilibration purposes
	inputFile >> LY >> endOfLine;				if (rank==ROOT) cout << "LY " << LY << endl;
	inputFile >> LZ >> endOfLine;				if (rank==ROOT) cout << "LZ " << LZ << endl;
	inputFile >> nbEqStep >> endOfLine;			if (rank==ROOT) cout << "nbEqStep " << nbEqStep << endl;
	inputFile >> equilTime >> endOfLine;		if (rank==ROOT) cout << "equilTime " << equilTime << endl;
	inputFile >> infoStep >> endOfLine;			if (rank==ROOT) cout << "infoStep " << infoStep <<endl;
	inputFile >> writeStep >> endOfLine;		if (rank==ROOT) cout << "writeStep " << writeStep << endl;
	inputFile >> floorTime >> endOfLine;		if (rank==ROOT) cout << "floorTime " << floorTime << endl;
	inputFile >> tetac >> endOfLine;			if (rank==ROOT) cout << "teta1  " << tetac  << endl;
	inputFile >> teta2 >> endOfLine;			if (rank==ROOT) cout << "teta2  " << teta2  << endl;
	inputFile >> teta_CB_fix >> endOfLine;		if (rank==ROOT) cout << "teta_CB_fix  " << teta_CB_fix  << endl;
	inputFile >> tauliquid >> endOfLine;		if (rank==ROOT) cout << "tauliquid " << tauliquid << endl;
	inputFile >> taugas >> endOfLine;			if (rank==ROOT) cout << "taugas " << taugas << endl;
	inputFile >> tau1  >> endOfLine;			if (rank==ROOT) cout << "tau_p " << tau1 << endl;
	inputFile >> tausurf >> endOfLine;	        if (rank==ROOT) cout << "tausurf " << tausurf << endl;
	inputFile >> kappa_p >> endOfLine;			if (rank==ROOT) cout << "kappa_p " << kappa_p << endl;
	inputFile >> A >> B >> gama >> endOfLine;	if (rank==ROOT) cout << "A, B, gama " << A << " " << B << " " << gama << endl;
	inputFile >> G[0] >> G[1] >> G[2] >> endOfLine; if (rank==ROOT) cout << "G " << G[0] << " " << G[1] << " " << G[2] << endl;
	inputFile >> p_thresh  >> endOfLine;		if (rank==ROOT) cout << "p_thresh " << p_thresh << endl;
	inputFile >> folder >> endOfLine;           if (rank==ROOT) cout << "Folder name is " << folder << endl;
    
    
    
    inputFile.close();
	tetac = tetac/180.0*M_PI; // Conversion from deg to rad
	teta2 = teta2/180.0*M_PI; // Conversion from deg to rad
	//cout << "teta_CB_fix deg" << teta_CB_fix << endl;
	teta_CB_fix = teta_CB_fix/180.0*M_PI; // Conversion from deg to rad
	//cout << "teta_CB_fix rad" << teta_CB_fix << endl;
	
	ifstream inputFile2("substrate.par");
	if (!inputFile2) {
		cout << "Can't open the input file substrate.par " << endl;
	}
	inputFile2 >> geometry >> endOfLine;		if (rank==ROOT) cout << "geometry " << geometry << endl;
	inputFile2 >> Dx >> endOfLine;				if (rank==ROOT) cout << "Dx " << Dx << endl;
	inputFile2 >> Dy >> endOfLine;				if (rank==ROOT) cout << "Dy " << Dy << endl;
	inputFile2 >> Dh >> endOfLine;				if (rank==ROOT) cout << "Dh " << Dh << endl;
	inputFile2 >> PeriodX >> endOfLine;			if (rank==ROOT) cout << "PeriodX " << PeriodX << endl;
	inputFile2 >> PeriodY >> endOfLine;			if (rank==ROOT) cout << "PeriodY " << PeriodY << endl;
	inputFile2 >> PostSymm >> endOfLine;		if (rank==ROOT) cout << "PostSymm " << PostSymm << endl;
    inputFile2 >> Px >> endOfLine;		if (rank==ROOT) cout << "Px " << Px << endl;
    inputFile2 >> Py >> endOfLine;		if (rank==ROOT) cout << "Py " << Py << endl;

	//inputFile2 >> gapWidth >> endOfLine; if (rank==ROOT) cout << "gapWidth " << gapWidth << endl;
	inputFile2.close();	 
	
	
	
	ifstream inputFile3("spherical.par");
	if (!inputFile3) {
		cout << "Can't open the input file spherical.par" << endl;
	}
	inputFile3 >> dropletR >> endOfLine;		if (rank==ROOT) cout << "dropletR " << dropletR << endl;// use -1 for no droplet
	inputFile3 >> dropletCenterX >> endOfLine;	if (rank==ROOT) cout << "dropletCenterX " << dropletCenterX << endl;
	inputFile3 >> dropletCenterY >> endOfLine;	if (rank==ROOT) cout << "dropletCenterY " << dropletCenterY << endl; // always use 1 in 3D simulation using symmetry
	inputFile3 >> dropletCenterZ >> endOfLine;	if (rank==ROOT) cout << "dropletCenterZ " << dropletCenterZ << endl;
	inputFile3 >> initUX >> initUY >> initUZ >> endOfLine; if (rank==ROOT) cout << "initU x,y,z, " << initUX << " " << initUY << " " << initUZ << endl;
	inputFile3.close();
	 
}

/*int wet::ReadInput(void)
 {
 
 string endOfLine;				//To store the one-line-comment in wet.dat 
	
	cout << "Process " << rank << ": reading input..." << endl;

	ifstream inputFile("wet.dat");
	if (!inputFile) 
	{
		cout << "Can't open the input file wet.dat" << endl;
		return 1;
	}
	inputFile >> LX >> endOfLine;  			//Number of X lattice points, read in wet.dat
	inputFile >> LY >> endOfLine;   		//Number of Y lattice points, read in wet.dat
	inputFile >> LZ >> endOfLine;  			//Number of Z lattice points, read in wet.dat
	inputFile >> nbEqStep >> endOfLine;   		//Number of total iterations, read in wet.dat
	inputFile >> infoStep >> endOfLine;   		//Number of steps at witch info is collected, read in wet.dat
	inputFile >> writeStep >> endOfLine;   		//Number of steps at which data is written, read in wet.dat
	inputFile >> teta1 >> endOfLine;  		//Planar surface contact angles of lower wall, read in wet.dat 
	inputFile >> teta2 >> endOfLine;  		//Planar surface contact angles of upper walls, read in wet.dat 
	inputFile >> tauliquid >> endOfLine;  		//Relaxation times of liquid, connected with viscosity (see literature), read in wet.dat
	inputFile >> taugas >> endOfLine;  		//Relaxation times of gas, connected with viscosity (see literature), read in wet.dat
	inputFile >> tau1  >> endOfLine;  		//Relaxation times of gas, connected with viscosity (see literature), read in wet.dat
	inputFile >> kappa_p >> endOfLine; 		//Proportional costant of distorsion term in free energy, (see literature), read in wet.dat
	inputFile >> A >> gama >> endOfLine;  		//Proportional costant of expansion term in free energy, and diffusion, read in wet.dat
	inputFile >> G[0] >> G[1] >> G[2] >> endOfLine;	//Force on the fluid
	inputFile.close();

	teta1=teta1/180.0*M_PI; // Conversion from deg to rad
	teta2=teta2/180.0*M_PI; // Conversion from deg to rad


	ifstream inputFile2("rectposts.dat");
	if (!inputFile2) 
	{
		cout << "Can't open the input file rectposts.dat (needed for floor and gap data as well, so always read)" << endl;
	}
	inputFile2 >> Dx >> endOfLine;  		//X dimension of a reactangular post (in lattice points), used in initialiseSurface();
	inputFile2 >> Dy >> endOfLine; 			//Y dimension of a reactangular post (in lattice points), used in initialiseSurface();
	inputFile2 >> Dh >> endOfLine; 			//Z dimension of a reactangular post (in lattice points), used in initialiseSurface();
	inputFile2 >> Period >> endOfLine;  		//Period of posts (equal for x and y dimension), used in initialiseSurface();
	inputFile2 >> PostSymm >> endOfLine;  		//Define offset of post beginning from x=0, used in initialiseSurface();
	inputFile2.close();	 

	
	
	ifstream inputFile3("spherical.dat");
	if (!inputFile3) 
	{
		cout << "Can't open the input file spherical.dat" << endl;
	}
	inputFile3 >> dropletR >> endOfLine;   		//Radius of a free droplet, used to calculate the new geometry of deformed drop
	inputFile3 >> dropletCenterX >> endOfLine;  	//X position of free droplet (in lattice points), also for deformated drop
	inputFile3 >> dropletCenterY >> endOfLine;  	//Y position of free droplet (in lattice points), also for deformated drop
	inputFile3 >> dropletCenterZ >> endOfLine;  	//Z position of free droplet (in lattice points), defined in spherical.dat
	inputFile3 >> initUX >> initUY >> initUZ >> endOfLine;  //Initial velocity of the droplet, used in LGConfig()
	inputFile3.close();	

	cout << "Process " << rank << ": input read." << endl;

} 
*/
