//writeJumpFile.cpp: description
//Lisa, 12 Aug 2011

#include "wet.h"

void wet::writeJumpFile(void)
{	
	char filenameic[10];
    string filenamei;
    snprintf(filenameic,10, "/jump.dat");
	filenamei=folder+filenameic;
	ofstream file6(filenamei.c_str(), ios::app);
	file6.precision(12);
	if (rank == ROOT)        file6 << t << "   "<< teta1*180/M_PI  << " " << teta_CB*180/M_PI  << " " << areaFrac << "   " << kin <<  " "<< kin_CM << " "<< RNtot << "     " << surfArea_init << " "                    << clLength_n_init << "     " << energy_init << " " << bulkE_init << " "<< interfaceE_init << " "<< surfaceE_init << "    " << energy_n_init << " " << bulkE_n_init << " "<< interfaceE_n_init << " "<< surfaceE_n_init << "        " << energy<< " " << bulkE<< " "<< interfaceE<< " "<< surfaceE<< "    " << energy_n<< " " << bulkE_n << " "<< interfaceE_n << " "<< surfaceE_n << "    " << eccentricity << endl ;
	//old: if (rank == ROOT) file6 << t << "   "<< teta1*180/M_PI  << " " << teta_CB*180/M_PI  << " " << areaFrac << "   " << kin <<  " "<< kin_CM <<                "     " << surfArea_init << " " << clLength << " " << clLength_n_init << "     " << energy_init << " " << bulkE_init << " "<< interfaceE_init << " "<< surfaceE_init << "    " << energy_n_init << " " << bulkE_n_init << " "<< interfaceE_n_init << " "<< surfaceE_n_init << "        " << energy<< " " << bulkE<< " "<< interfaceE<< " "<< surfaceE<< "    " << energy_n<< " " << bulkE_n << " "<< interfaceE_n << " "<< surfaceE_n << "    " << eccentricity << endl;
	
	
	file6.close();
	
	
    char filenameic2[10];
    snprintf(filenameic2,10, "/jump2.dat");
	filenamei=folder+filenameic;
	ofstream file7(filenamei.c_str(), ios::app);
	file7.precision(12);  
	//if (rank == ROOT) file7 << equilfirst << " " << duplicationtype << " " << dimensions << " " << LX << " " << LY << " " << LZ << " " << equilTime << " "                     << teta1*180/M_PI  << " " << tauliquid << " " << taugas << " " << kappa_p << " " << A << " " << gama << " " << G[0] << " " << G[1] << " " << G[2] << " " << p_thresh << "   " << geometry << " " << Dx << " " << Dy << " " << Dh << " " << PeriodX << " " << PeriodY << " " << PostSymm << "   " << dropletR << " " << dropletCenterX << " " << dropletCenterY << " " << dropletCenterZ << " " << initUX << " " << initUY << " " << initUZ << "       " << equilTime << " " << t << " " << teta_CB*180/M_PI << " " << teta_CB_fix*180/M_PI  << " " << areaFrac << " "                                  << kinL << " " << kin_CM << " " << RVZcm << " "  << RNtot << " " << eccentricity << " " << surfArea_init << " " << clLength_n_init << "     " << energy_n_init << " " << bulkE_n_init << " "<< interfaceE_n_init << " "<< surfaceE_n_init << " " << energy_n << " " << bulkE_n << " "<< interfaceE_n << " "<< surfaceE_n  << endl;
	  if (rank == ROOT) file7 << equilfirst << " " << duplicationtype << " " << dimensions << " " << LX << " " << LY << " " << LZ << " " << equilTime << " " << floorTime << " " << teta1*180/M_PI  << " " << tauliquid << " " << taugas << " " << kappa_p << " " << A << " " << gama << " " << G[0] << " " << G[1] << " " << G[2] << " " << p_thresh << "   " << geometry << " " << Dx << " " << Dy << " " << Dh << " " << PeriodX << " " << PeriodY << " " << PostSymm << "   " << dropletR << " " << dropletCenterX << " " << dropletCenterY << " " << dropletCenterZ << " " << initUX << " " << initUY << " " << initUZ << "       " << equilTime << " " << t << " " << teta_CB*180/M_PI << " " << teta_CB_fix*180/M_PI  << " " << areaFrac << " " << kin_init << " " << kin << " " << kinL << " " << kin_CM << " " << RVZcm << " "  << RNtot << " " << eccentricity << " " << surfArea_init << " " << clLength_n_init << "     " << energy_n_init << " " << bulkE_n_init << " "<< interfaceE_n_init << " "<< surfaceE_n_init << " " << energy_n << " " << bulkE_n << " "<< interfaceE_n << " "<< surfaceE_n  << "         " << dissInt_liq << " " << dissInt_int << " " << dissInt_gas << " " << diffDissInt_liq << " " << diffDissInt_int<< " " << diffDissInt_gas << " " << mLiq << " " << mInt << " " << mGas << endl;

	file7.close();

}

