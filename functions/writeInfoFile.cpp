//saveFiles.cpp: save data in files
//Lisa, 12 Aug 2011

#include "wet.h"



void wet::writeInfoFile(void)
{
    char filenamei[50];
    snprintf(filenamei,50, "./uxs%.3frat%.3foff%.0f/info.dat",initUX,dropletR/Dx,dropletCenterY-Px-Dx/2);
	String fileName5(filenamei);
	ofstream file5(fileName5.get(), ios::app);
	
	if (t==0 && rank==ROOT){
		file5.precision(4);
		
		
		file5 << "Simulation settings " << endl;
		file5 << " " << endl;
		file5 << "equilfirst? " << equilfirst << endl;
		file5 << "velocityprofile? " << velocityprofile << endl;
		file5 << "duplicationtype " << duplicationtype  << endl;
		file5 << " " << endl;
		
		file5 << "Simulation parameters " << endl;
		file5 << " " << endl;
		file5 << "dimensions " << dimensions << endl;
		file5 << "LX " << LXc << endl;   //Lisa: LXc is the constant value that is not to be changed. LX is changed for equilibration purposes
		file5 << "LY " << LY << endl;
		file5 << "LZ " << LZ << endl;
		file5 << "nbEqStep " << nbEqStep << endl;
		file5 << "equilTime " << equilTime << endl;
		file5 << "infoStep " << infoStep <<endl;
		file5 << "writeStep " << writeStep << endl;
		file5 << "floorTime " << floorTime << endl;
		file5 << "teta1  " << tetac  << " = " << tetac*180/M_PI << endl;
		file5 << "teta2  " << teta2  << " = " << teta2*180/M_PI << endl;
		file5 << "teta_CB_fix  " << teta_CB_fix  << " = " << teta_CB_fix*180/M_PI << endl;
		file5 << "tauliquid " << tauliquid << endl;
		file5 << "taugas " << taugas << endl;
		file5 << "tau_p " << tau1 << endl;
		file5 << "tausurf " << tausurf << endl;
		file5 << "kappa_p " << kappa_p << endl;
		file5 << "A, B, gama " << A << " " << B << " " << gama << endl;
		file5 << "G " << G[0] << " " << G[1] << " " << G[2] << endl;
		file5 << "p_thresh " << p_thresh <<" p_thresh_natural " << tanh(1.0)<< endl;
		file5 << "  " << endl;
		
		file5 << "Solid parameters " << endl;
		file5 << " " << endl;
		file5 << "geometry " << geometry << endl;
		file5 << "Dx " << Dx << endl;
		file5 << "Dy " << Dy << endl;
		file5 << "Dh " << Dh << endl;
		file5 << "PeriodX " << PeriodX  << endl;
		file5 << "PeriodY " << PeriodY  << endl;
		file5 << "PostSymm " << PostSymm  << endl;
		file5 << "Px "  << Px << endl;
        file5 << "Py "  <<Py << endl;
		file5 << "gapWidth "  <<endl;
        file5 << "  " << endl;
        
        file5 << "Calculated Experimental Touch parameters" << endl;
        file5 << "Surface tension Sl " << stl*mo/(dt*dt) << endl;
        file5 << "Surafce tension lg " << sqrt(8.0*kappa_p/A)*mo/(dt*dt) << endl;
        file5 << "Surface Tension gas solid " << stg*mo/(dt*dt) << endl;
        file5 << "Liquid Viscosity " << vl*c*c*dt << endl;
        file5 << "Gas viscosity " << vg*c*c*dt << endl;
        file5 << "Reynolds number " << Re << endl;
        file5 << "Weber number " << we << endl;
        file5 << "Oh number " << Oh << endl;
        file5 << "wetting potential " << w << endl << endl;
		
		file5 << "Drop parameters " << endl;
		file5 << " " << endl;
		file5 << "dropletR " << dropletR << endl;
		file5 << "dropletCenterX " << dropletCenterX  << endl;
		file5 << "dropletCenterY " << dropletCenterY  << endl;
		file5 << "dropletCenterZ " << dropletCenterZ  << endl;
		file5 << "initU x,y,z, " << initUX << initUY<< initUZ << endl;
		file5 << "  " << endl;
		file5 << " ---------------------------------------- " << endl;
		
		file5 << "Calculated quantities" << endl;
		file5 << " " << endl;
		file5 << "Interface width = sqrt(8*kappa/A) " << sqrt(8.0*kappa_p/A) << endl;
		file5 << "Surface tension = sqrt(8/9*kappa*A) " << sqrt(8.0/9.0*kappa_p*A) << endl;
		file5 << "CB angle " << teta_CB*180.0/M_PI << endl;
		if (teta_CB_fix != 0) {
			file5 << "which was fixed externally " << endl;
		}
		file5 << "Area fraction " << areaFrac << endl;
		file5 << " " << endl;
	}
	
	if (t==equilTime-1 && rank==ROOT){
		file5.precision(12);
		file5 << "BEFORE drop duplication, t = " << t << endl;
		file5 << "totalF      ( .., _n) " << energy  << "  " << energy_n<< endl;
		file5 << "bulkF       ( .., _n) " << bulkE   << "  " << bulkE_n << endl;
		file5 << "surfaceF    ( .., _n) " << surfaceE << "  " << surfaceE_n << endl;
		file5 << "interfaceF  ( .., _n) " << interfaceE << "  " << interfaceE_n<< endl;
		file5 << "length of CL ( .., _n) "<< clLength << "  " << clLength_n << endl;
		file5 << "contact area " << surfArea << endl; 		file5 << " " << endl;
	}
	if (t==equilTime && rank==ROOT){
		file5.precision(12);
		//file5 << "AFTER drop duplication, t = " << t << endl;
		file5 << "totalF      ( .., _n) " << energy  << "  " << energy_n<< endl;
		file5 << "bulkF       ( .., _n) " << bulkE   << "  " << bulkE_n << endl;
		file5 << "surfaceF    ( .., _n) " << surfaceE << "  " << surfaceE_n << endl;
		file5 << "interfaceF  ( .., _n) " << interfaceE << "  " << interfaceE_n<< endl;
		file5 << "length of CL ( .., _n) "<< clLength << "  " << clLength_n << endl;
		file5 << "contact area " << surfArea << endl;
		file5 << " " << endl;
	}
	
	/*if (t>equilTime && detachflag == 1 && rank==ROOT){
     file5.precision(12);
     file5 << "--- DROP DETACHED! --- t  = " << t << endl;
     file5 << "totalF      ( .., _n) " << energy  << "  " << energy_n<< endl;
     file5 << "bulkF       ( .., _n) " << bulkE   << "  " << bulkE_n << endl;
     file5 << "surfaceF    ( .., _n) " << surfaceE << "  " << surfaceE_n << endl;
     file5 << "interfaceF  ( .., _n) " << interfaceE << "  " << interfaceE_n<< endl;
     file5 << "length of CL ( .., _n) "<< clLength << "  " << clLength_n << endl;
     file5 << "contact area " << surfArea << endl;
     }*/
	file5.close();	
}
