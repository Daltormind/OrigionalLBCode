//ARolling.cpp: this function computes dynamic variables (kinetic energy, centre of mass velocity,...)
//Marco, 29 June 2011, Lisa, 26 Sept 2011, changed all weighting to n[k] instead of p[k] 
//Lisa, 01 Oct 2011, included factor 1/2 in kinetic energy (Ekin = 1/2(!!!)mv^2)

#include "wet.h"

void wet::ARolling(void)
{  
	//cout << " *** A Rolling started" << endl;
	
	double ECCX, ECCZ; //eccentricity of drop in x and y direction
	double vx_sq, vz_sq;
	double vmax = 0.0;
	

	

	ECCX = 0.0L;
	ECCZ = 0.0L;	

	RNtot = 0.0L;
	RNtot_gas = 0.0L;
	
	RXcm=0.0L; RZcm=0.0L;
	RVXcm=0.0L; RVYcm=0.0L; RVZcm=0.0L;
	
	kinL = 0.0L; //kinetic energy of liq (incl interface)
	kinG = 0.0L;
	kin = 0.0L;
	
	vx_sq = 0.0L;
	vz_sq = 0.0L;

	for(k = k1; k < k2; k++)
	{
		computeCoordinate(k);
		
		d_xk = double(xk);
		d_yk = double(yk);
		d_zk = double(zk);
		
		//if (mask[k] != 1){   // only for testing exclusion of some areas in dissipation calculation (ARolling, computeDissipation, computeFreeEnergy)
		//computeCoordinate(k);
		//if (zk > Dh + dropletCenterZ && zk + dropletCenterZ < LZ){	
			
			if (p[k]> 0.0L && mask[k]!=28 ) // if we are in the droplet
			{					
				RNtot += n[k]; 			//sum the amount of mass, is it better phi?
				RXcm  += n[k]*d_xk; 		//calculation of x coordinate of cantre of mass of fluid f
				RZcm  += n[k]*d_zk;		//calculation of z coordinate of cantre of mass of fluid f
				
				RVXcm += n[k]*uxs[k]; 	//calculation of velocity of center of mass in x direction (remember to compare RVX d RX/dt) 
				RVYcm += n[k]*uys[k]; 	//calculation of velocity of center of mass in y direction, not relevant in 2d
				RVZcm += n[k]*uzs[k];	//calculation of velocity of center of mass in z direction
				
				kinL += 0.5*n[k]*(uxs[k]*uxs[k] + uys[k]*uys[k] + uzs[k]*uzs[k]);	//TOTAL kinetic energy of drop
				
				vx_sq += uxs[k]*uxs[k];
				vz_sq += uzs[k]*uzs[k];
				
			} 
			if (p[k]<= 0.0L && mask[k]!=28 ) // outside the droplet, aka in the gas
			{
				RNtot_gas += n[k];
				kinG += 0.5*n[k]*(uxs[k]*uxs[k] + uys[k]*uys[k] + uzs[k]*uzs[k]);
				
				if (uxs[k]*uxs[k] + uys[k]*uys[k] + uzs[k]*uzs[k] > vmax)
				{
					vmax = sqrt(uxs[k]*uxs[k] + uys[k]*uys[k] + uzs[k]*uzs[k]);
				}
			}
			
		//}

	} 
	
	

	// Summing all the variables and putting that on ROOT
	
	double tempvar = 0.0;
	MPI_Reduce(&RNtot,&tempvar,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	RNtot = tempvar;
	MPI_Reduce(&RNtot_gas,&tempvar,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	RNtot_gas = tempvar;

	MPI_Reduce(&RXcm,&tempvar,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	RXcm = tempvar;
	MPI_Reduce(&RZcm,&tempvar,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	RZcm = tempvar;

	MPI_Reduce(&RVXcm,&tempvar,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	RVXcm = tempvar;
	MPI_Reduce(&RVYcm,&tempvar,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	RVYcm = tempvar;
	MPI_Reduce(&RVZcm,&tempvar,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	RVZcm = tempvar;
	
	MPI_Reduce(&kinL,&tempvar,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	kinL = tempvar;
	MPI_Reduce(&kinG,&tempvar,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	kinG = tempvar;
	
	MPI_Reduce(&vx_sq,&tempvar,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	vx_sq = tempvar;
	MPI_Reduce(&vz_sq,&tempvar,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	vz_sq = tempvar;
	MPI_Reduce(&RNtot_gas,&tempvar,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	RNtot_gas = tempvar;
	
	
	MPI_Reduce(&vmax,&tempvar,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD); //max of all processes' vmax
	vmax = tempvar;
	
		
	if (rank == ROOT) {
		
		RXcm /= RNtot; 
		RZcm /= RNtot;
		
		kin_CM = 0.5*( RVXcm*RVXcm + RVYcm*RVYcm + RVZcm*RVZcm );
		kin_CM /= RNtot;
		
		ratio_kin = kin_CM / kinL ; 
		
		kin = kinL + kinG;
		
		RVXcm /= RNtot;
		RVZcm /= RNtot;  
		
		vx_sq /= RNtot;
		vz_sq /= RNtot;
		
	}
	
	
	//  calculate eccentricity: distribute centre of mass coordinate back to all processes, calc eccentricity, put it on ROOT
	
	MPI_Bcast(&RXcm,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&RZcm,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	
	//cout << "Process " << rank << " : RXcm = " << RXcm << " RZcm = " << RZcm << endl;
	
	for (k = k1; k < k2; k++){
		
		if (p[k]> 0.0L && mask[k]!=28 ) {
			computeCoordinate(k);
			if (xk == int(RXcm)) ECCZ += n[k];
			if (zk == int(RZcm)) ECCX += n[k];
		}
	}
	
	MPI_Reduce(&ECCX,&tempvar,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	ECCX = tempvar;
	MPI_Reduce(&ECCZ,&tempvar,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	ECCZ = tempvar;
	
	//cout << "t, ECCX " << t << " " << ECCX << endl;
	//cout << "t, ECCZ " << t << " " << ECCZ << endl;
	
	
	 
	if (rank == ROOT) {
		
		eccentricity = ECCX/ECCZ;
		
		char filenameic[20];
        string filenamei;
        snprintf(filenameic,20, "/totalmomenta.dat");
        filenamei=folder+filenameic;
		ofstream file3(filenamei.c_str(), ios::app);
		file3.precision(7);
		//file3<< t <<" " << RXcm<<" "<<RVXcm<<" "<<RXcm1+0.4244132*dropletR*sqrt2 <<" "<<RVXcm1<<" "<<RVZcm<<" "<< ECCX << " "<<kin<<" "<<" "<<kin_CM<<endl; 
		file3<< t <<" "<< kin_CM << " "<< surfArea <<"     " << RXcm<<" "<<RVXcm<<" "<<RZcm<<" "<< RVZcm << "     " << eccentricity <<" "<<kinL<<" "<<kin_CM<< " " << RNtot  << "      " << vx_sq << " " << vz_sq <<"      " << kinG << " " << RNtot_gas << "   " << vmax << endl; 
		file3.close();
		
	}
		
}

