//ReadInput.cpp: read input from wet.dat, rectpost.dat and spherical.dat
//Marco, 28 June 2011, Lisa, 28 July 2011

#include "wet.h"


int wet::ReadInput(void)
{	
	
	 
	ifstream inputFile;
//---------------------------Get inputs from the settings file---------------
    inputFile.open("settings.par");
	
    if (!inputFile) {
		cout << "Can't open the input file settings.par " << endl;
	}
	
    inputFile >> equilfirst ;
    inputFile.ignore(250,'\n');
    if (rank==ROOT) cout << "equilfirst? " << equilfirst << endl;
	
    inputFile >> velocityprofile ;
    inputFile.ignore(250,'\n');
    if (rank==ROOT) cout << "velocityprofile? " << velocityprofile << endl;
	
    inputFile.close();

//---------------------- Get inputs from wet.par ----------------------------
    
    inputFile.open("wet.par");
	if (!inputFile) {
		cout << "Can't open the input file wet.par" << endl;
		return 1;
	}
	
	inputFile >> dimensions;
    inputFile.ignore(250,'\n');
    if (rank==ROOT) cout << "dimensions " << dimensions << endl;
	if (dimensions != 2 && dimensions !=3) {
		cout << "Dimensions not well defined in settings.par (2 or 3?)" << endl;
		return 1;
	}
	
    inputFile >> LX ;
    inputFile.ignore(250,'\n');
    if (rank==ROOT) cout << "LX " << LX << endl;   
	
    inputFile >> LY ;
    inputFile.ignore(250,'\n');
    if (rank==ROOT) cout << "LY " << LY << endl;
	
    inputFile >> LZ ;
    inputFile.ignore(250,'\n');
    if (rank==ROOT) cout << "LZ " << LZ << endl;
	
    inputFile >> nbEqStep ;
    inputFile.ignore(250,'\n');
    if (rank==ROOT) cout << "nbEqStep " << nbEqStep << endl;
	
    inputFile >> equilTime ;
    inputFile.ignore(250,'\n');
    if (rank==ROOT) cout << "equilTime " << equilTime << endl;
	
    inputFile >> infoStep ;
    inputFile.ignore(250,'\n');
    if (rank==ROOT) cout << "infoStep " << infoStep <<endl;
	
    inputFile >> writeStep ;
    inputFile.ignore(250,'\n');
    if (rank==ROOT) cout << "writeStep " << writeStep << endl;
	
    inputFile >> floorTime ;
    inputFile.ignore(250,'\n');
    if (rank==ROOT) cout << "floorTime " << floorTime << endl;
	
    inputFile >> tetac ;
    inputFile.ignore(250,'\n');
    if (rank==ROOT) cout << "teta1  " << tetac  << endl;
	
    inputFile >> teta2 ;
    inputFile.ignore(250,'\n');
    if (rank==ROOT) cout << "teta2  " << teta2  << endl;
	
    inputFile >> teta_CB_fix ;
    inputFile.ignore(250,'\n');
    if (rank==ROOT) cout << "teta_CB_fix  " << teta_CB_fix  << endl;
	
    inputFile >> tauliquid ;
    inputFile.ignore(250,'\n');
    if (rank==ROOT) cout << "tauliquid " << tauliquid << endl;
	
    inputFile >> taugas;
    inputFile.ignore(250,'\n');
    if (rank==ROOT) cout << "taugas " << taugas << endl;
	
    inputFile >> tau1  ;
    inputFile.ignore(250,'\n');
    if (rank==ROOT) cout << "tau_p " << tau1 << endl;
	
    inputFile >> tausurf ;
    inputFile.ignore(250,'\n');
    if (rank==ROOT) cout << "tausurf " << tausurf << endl;
	
    inputFile >> kappa_p;
    inputFile.ignore(250,'\n');
    if (rank==ROOT) cout << "kappa_p " << kappa_p << endl;
	
    inputFile >> A >> B >> gama ;
    inputFile.ignore(250,'\n');
    if (rank==ROOT) cout << "A, B, gama " << A << " " << B << " " << gama << endl;
	
    inputFile >> G[0] >> G[1] >> G[2] ;
    inputFile.ignore(250,'\n');
    if (rank==ROOT) cout << "G " << G[0] << " " << G[1] << " " << G[2] << endl;
	
    inputFile >> p_thresh  ;
    inputFile.ignore(250,'\n');
    if (rank==ROOT) cout << "p_thresh " << p_thresh << endl;
	
    inputFile >> folder ;
    inputFile.ignore(250,'\n');
    if (rank==ROOT) cout << "Folder name is " << folder << endl;
    
    
    
    inputFile.close();
	tetac = tetac/180.0*M_PI; // Conversion from deg to rad
	teta2 = teta2/180.0*M_PI; // Conversion from deg to rad
	
	teta_CB_fix = teta_CB_fix/180.0*M_PI; // Conversion from deg to rad

	//-------------------Inputs from substrate.par------------------------------
	
    inputFile.open("substrate.par");
	
    if (!inputFile) {
		cout << "Can't open the input file substrate.par " << endl;
	}
	
    inputFile >> geometry ;
    inputFile.ignore(250,'\n');
    if (rank==ROOT) cout << "geometry " << geometry << endl;
	
    inputFile >> Dx ;
    inputFile.ignore(250,'\n');
    if (rank==ROOT) cout << "Dx " << Dx << endl;
	
    inputFile >> Dy ;
    inputFile.ignore(250,'\n');
    if (rank==ROOT) cout << "Dy " << Dy << endl;
	
    inputFile >> Dh ;
    inputFile.ignore(250,'\n');
    if (rank==ROOT) cout << "Dh " << Dh << endl;
	
    inputFile >> PeriodX ;
    inputFile.ignore(250,'\n');
    if (rank==ROOT) cout << "PeriodX " << PeriodX << endl;

    inputFile >> PeriodY ;
    inputFile.ignore(250,'\n');
    if (rank==ROOT) cout << "PeriodY " << PeriodY << endl;
	
    inputFile >> PostSymm ;
    inputFile.ignore(250,'\n');
    if (rank==ROOT) cout << "PostSymm " << PostSymm << endl;
    
    inputFile >> Px ;
    inputFile.ignore(250,'\n');
    if (rank==ROOT) cout << "Px " << Px << endl;
    
    inputFile >> Py ;
    inputFile.ignore(250,'\n');
    if (rank==ROOT) cout << "Py " << Py << endl;
    
    inputFile >> Pz ;
    inputFile.ignore(250,'\n');
    if (rank==ROOT) cout << "Pz " << Pz << endl;

	
	inputFile.close();	 
	
//-------------------------Input from spherical.par--------------------------------
	
    inputFile.open("spherical.par");
	if (!inputFile) {
		cout << "Can't open the input file spherical.par" << endl;
	}
	
    inputFile >> dropletR ;
    inputFile.ignore(250,'\n');
    if (rank==ROOT) cout << "dropletR " << dropletR << endl;// use -1 for no droplet
	
    inputFile >> dropletCenterX ;
    inputFile.ignore(250,'\n');
    if (rank==ROOT) cout << "dropletCenterX " << dropletCenterX << endl;
	
    inputFile >> dropletCenterY;
    inputFile.ignore(250,'\n');
    if (rank==ROOT) cout << "dropletCenterY " << dropletCenterY << endl; // always use 1 in 3D simulation using symmetry
	
    inputFile >> dropletCenterZ ;
    inputFile.ignore(250,'\n');
    if (rank==ROOT) cout << "dropletCenterZ " << dropletCenterZ << endl;
	
    inputFile >> initUX >> initUY >> initUZ ;
    inputFile.ignore(250,'\n');
    if (rank==ROOT) cout << "initU x,y,z, " << initUX << " " << initUY << " " << initUZ << endl;
	
    inputFile.close();
	 
    
    return(1);
}

