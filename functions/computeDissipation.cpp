
#include "wet.h"



void wet::computeDissipation(void)
{	
	double dux_dx,dux_dy,dux_dz;
	double duy_dx,duy_dy,duy_dz;
	double duz_dx,duz_dy,duz_dz;
	double dmu_dx, dmu_dy, dmu_dz;
	
	double energyDiff; 
	double kinDiff;
		
	double divV;
	double kin_liq, kin_int, kin_gas;
	
	double tau22;
	
	double diffDiss_surf, diss_surf;
	
	//double checkEBalDiff = 0.0;
	
	
	dissSum_liq = 0.0;
	dissSum_int = 0.0;
	dissSum_gas = 0.0;
	diffDissSum_liq = 0.0;
	diffDissSum_int = 0.0;
	diffDissSum_gas = 0.0;
		
	mLiq = 0.0;
	mGas = 0.0;
	mInt = 0.0;
	
	kin_liq = 0.0L;
	kin_int = 0.0L;
	kin_gas = 0.0L;
	
	
	diffDiss_surf = 0.0L;
	diss_surf = 0.0L;
	
	//cout << "Process " << rank << ": computeDiss check1, t = "<< t << endl;
	
	for (k=k1; k<k2; k++) {
	
		if (mask[k] != 28){
			
			dux_dx = ((uxs[dd1[k]] - uxs[dd2[k]])/6 + (uxs[dda[k]]+uxs[ddc[k]]+uxs[ddi[k]]+uxs[ddk[k]]-uxs[ddb[k]]-uxs[ddd[k]]-uxs[ddj[k]]-uxs[ddl[k]])/12);
			dux_dy = ((uxs[dd3[k]] - uxs[dd4[k]])/6 + (uxs[dda[k]]+uxs[ddb[k]]+uxs[dde[k]]+uxs[ddg[k]]-uxs[ddc[k]]-uxs[ddd[k]]-uxs[ddf[k]]-uxs[ddh[k]])/12);
			dux_dz = ((uxs[dd5[k]] - uxs[dd6[k]])/6 + (uxs[dde[k]]+uxs[ddf[k]]+uxs[ddi[k]]+uxs[ddj[k]]-uxs[ddg[k]]-uxs[ddh[k]]-uxs[ddk[k]]-uxs[ddl[k]])/12);
			
			duy_dx = ((uys[dd1[k]] - uys[dd2[k]])/6 + (uys[dda[k]]+uys[ddc[k]]+uys[ddi[k]]+uys[ddk[k]]-uys[ddb[k]]-uys[ddd[k]]-uys[ddj[k]]-uys[ddl[k]])/12);
			duy_dy = ((uys[dd3[k]] - uys[dd4[k]])/6 + (uys[dda[k]]+uys[ddb[k]]+uys[dde[k]]+uys[ddg[k]]-uys[ddc[k]]-uys[ddd[k]]-uys[ddf[k]]-uys[ddh[k]])/12);
			duy_dz = ((uys[dd5[k]] - uys[dd6[k]])/6 + (uys[dde[k]]+uys[ddf[k]]+uys[ddi[k]]+uys[ddj[k]]-uys[ddg[k]]-uys[ddh[k]]-uys[ddk[k]]-uys[ddl[k]])/12);
			
			duz_dx = ((uzs[dd1[k]] - uzs[dd2[k]])/6 + (uzs[dda[k]]+uzs[ddc[k]]+uzs[ddi[k]]+uzs[ddk[k]]-uzs[ddb[k]]-uzs[ddd[k]]-uzs[ddj[k]]-uzs[ddl[k]])/12);
			duz_dy = ((uzs[dd3[k]] - uzs[dd4[k]])/6 + (uzs[dda[k]]+uzs[ddb[k]]+uzs[dde[k]]+uzs[ddg[k]]-uzs[ddc[k]]-uzs[ddd[k]]-uzs[ddf[k]]-uzs[ddh[k]])/12);
			duz_dz = ((uzs[dd5[k]] - uzs[dd6[k]])/6 + (uzs[dde[k]]+uzs[ddf[k]]+uzs[ddi[k]]+uzs[ddj[k]]-uzs[ddg[k]]-uzs[ddh[k]]-uzs[ddk[k]]-uzs[ddl[k]])/12);
			
			// grad(mu) calculation for diffusion dissipation
			dmu_dx = ((mu[dd1[k]] - mu[dd2[k]])/6 + (mu[dda[k]]+mu[ddc[k]]+mu[ddi[k]]+mu[ddk[k]]-mu[ddb[k]]-mu[ddd[k]]-mu[ddj[k]]-mu[ddl[k]])/12);
			dmu_dy = ((mu[dd3[k]] - mu[dd4[k]])/6 + (mu[dda[k]]+mu[ddb[k]]+mu[dde[k]]+mu[ddg[k]]-mu[ddc[k]]-mu[ddd[k]]-mu[ddf[k]]-mu[ddh[k]])/12);
			dmu_dz = ((mu[dd5[k]] - mu[dd6[k]])/6 + (mu[dde[k]]+mu[ddf[k]]+mu[ddi[k]]+mu[ddj[k]]-mu[ddg[k]]-mu[ddh[k]]-mu[ddk[k]]-mu[ddl[k]])/12);
			
			
			
			tau22 = 0.5*((tauliquid+taugas)+p[k]*(tauliquid-taugas));
			if (tau22  < 0.5) {
				tau22 = 0.5;
			}
			
			dissipation[k] = (c2*dt/3.0*(tau22-0.5))* n[k] * ( (dux_dy + duy_dx)*(dux_dy + duy_dx) + (dux_dz + duz_dx)*(dux_dz + duz_dx) + (duy_dz + duz_dy)*(duy_dz + duz_dy) + 2.0*(dux_dx*dux_dx + duy_dy*duy_dy + duz_dz*duz_dz) );
			
			diffDiss[k] = dt*gama*(tau1-0.5)*( dmu_dx*dmu_dx + dmu_dy*dmu_dy + dmu_dz*dmu_dz  );
			
			/*eBal[k] += 0.5*n[k]*(uxs[k]*uxs[k] + uys[k]*uys[k] + uzs[k]*uzs[k]);   // kinetic energy part, free En part in computeFreeEnergy2.cpp, set to zero at start of computeFreeEnergy2.cpp
			// now eBal is sum of TotalFreeE[k] and kin[k]; then take difference to previous values = dissipated energy
			eBalDiff[k]  = -(eBal[k] - eBal_old[k])/infoStep;
			eBal_old[k] = eBal[k];
			
			checkEBalDiff += eBalDiff[k]; 
			*/
			//if (mask[k] != 1 ){ // only for testing exclusion of some areas in dissipation calculation (ARolling, computeDissipation, computeFreeEnergy)
			//computeCoordinate(k);
			//if (zk > Dh + dropletCenterZ && zk + dropletCenterZ < LZ){	
				
				if (fabs(p[k]) > p_thresh_n && p[k]> 0.0L)
				{
					dissSum_liq  += dissipation[k];
					diffDissSum_liq  += diffDiss[k];
					mLiq += n[k];
					kin_liq += 0.5*n[k]*(uxs[k]*uxs[k] + uys[k]*uys[k] + uzs[k]*uzs[k]);	
					
					if (mask[k] == 1) {
						diffDiss_surf += diffDiss[k];
						diss_surf += dissipation[k]; 
					}
					//divV = dux_dx + duy_dy + duz_dz;
					//dissSum2 += (c2*dt/3.0*(tau22-0.5)) * n[k] * ( (dux_dy + duy_dx)*(dux_dy + duy_dx) + (dux_dz + duz_dx)*(dux_dz + duz_dx) + (duy_dz + duz_dy)*(duy_dz + duz_dy) +  2.0*( (dux_dx - 1.0/3.0*divV)*(dux_dx - 1.0/3.0*divV) + (duy_dy - 1.0/3.0*divV)*(duy_dy - 1.0/3.0*divV) + (duz_dz - 1.0/3.0*divV)*(duz_dz - 1.0/3.0*divV) ) ) ;
				}
				if (fabs(p[k]) > p_thresh_n && p[k]< 0.0L)
				{
					dissSum_gas += dissipation[k];
					diffDissSum_gas  += diffDiss[k];
					mGas += n[k];
					kin_gas += 0.5*n[k]*(uxs[k]*uxs[k] + uys[k]*uys[k] + uzs[k]*uzs[k]);
					
					if (mask[k] == 1) {
						diffDiss_surf += diffDiss[k];
						diss_surf += dissipation[k]; 
					}
				}
				if (fabs(p[k]) <= p_thresh_n)
				{
					dissSum_int += dissipation[k];
					diffDissSum_int  += diffDiss[k];
					mInt += n[k];
					kin_int += 0.5*n[k]*(uxs[k]*uxs[k] + uys[k]*uys[k] + uzs[k]*uzs[k]);
					
					if (mask[k] == 1) {
						diffDiss_surf += diffDiss[k];
						diss_surf += dissipation[k]; 
					}
				}
				//cout << "Process " << rank << ": , t = " << t << " eta , " << " k "<< k << " " << (c2*dt/3.0*(tau22-0.5)) << endl;
			//}
			
			/*
			cout << "Process " << rank << ": , t = " << t << " dux_dx, " << " k "<< k << " " << dux_dx << endl;
			cout << "Process " << rank << ": , t = " << t << " dux_dy, " << " k "<< k << " " << dux_dy << endl;
			cout << "Process " << rank << ": , t = " << t << " dux_dz, " << " k "<< k << " " << dux_dz << endl;
			cout << "Process " << rank << ": , t = " << t << " duy_dx, " << " k "<< k << " " << duy_dx << endl;
			cout << "Process " << rank << ": , t = " << t << " duy_dy, " << " k "<< k << " " << duy_dy << endl;
			cout << "Process " << rank << ": , t = " << t << " duy_dz, " << " k "<< k << " " << duy_dz << endl;
			cout << "Process " << rank << ": , t = " << t << " duz_dx, " << " k "<< k << " " << duz_dx << endl;
			cout << "Process " << rank << ": , t = " << t << " duz_dy, " << " k "<< k << " " << duz_dy << endl;
			cout << "Process " << rank << ": , t = " << t << " duz_dz, " << " k "<< k << " " << duz_dz << endl;
			 */
			
		}		
		
	}
	//cout << "Process " << rank << ": computeDiss check2, t = "<< t << endl;
	
	double tempvar = 0.0;
	MPI_Reduce(&dissSum_liq,&tempvar,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	dissSum_liq = tempvar;
	MPI_Reduce(&dissSum_gas,&tempvar,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	dissSum_gas = tempvar;	
	MPI_Reduce(&dissSum_int,&tempvar,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	dissSum_int = tempvar;
	
	MPI_Reduce(&diffDissSum_liq,&tempvar,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	diffDissSum_liq = tempvar;
	MPI_Reduce(&diffDissSum_gas,&tempvar,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	diffDissSum_gas = tempvar;	
	MPI_Reduce(&diffDissSum_int,&tempvar,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	diffDissSum_int = tempvar;
	
	MPI_Reduce(&mLiq,&tempvar,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	mLiq = tempvar;
	MPI_Reduce(&mGas,&tempvar,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	mGas = tempvar;
	MPI_Reduce(&mInt,&tempvar,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	mInt = tempvar;
	
	MPI_Reduce(&kin_liq,&tempvar,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	kin_liq = tempvar;
	MPI_Reduce(&kin_gas,&tempvar,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	kin_gas = tempvar;	
	MPI_Reduce(&kin_int,&tempvar,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	kin_int = tempvar;
	
	MPI_Reduce(&diffDiss_surf,&tempvar,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	diffDiss_surf = tempvar;
	MPI_Reduce(&diss_surf,&tempvar,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	diss_surf = tempvar;	
	
	
	
	//cout << "Process " << rank << ": computeDiss check2, t = "<< t << endl;
	
	double reducedEnergy = 0.0;
	MPI_Reduce(&bulkE_n,&reducedEnergy,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	bulkE_n = reducedEnergy;
	
	reducedEnergy = 0.0;
	MPI_Reduce(&interfaceE_n,&reducedEnergy,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	interfaceE_n = reducedEnergy;
	
	reducedEnergy = 0.0;
	MPI_Reduce(&surfaceE_n,&reducedEnergy,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	surfaceE_n = reducedEnergy;
	
	energy_n = bulkE_n + interfaceE_n + surfaceE_n;
	
	//MPI_Reduce(&checkEBalDiff,&tempvar,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	//checkEBalDiff = tempvar;
	
	
	
	if (rank == ROOT){
		
		dissInt_liq = dissSum_liq*infoStep + dissInt_liq_old;
		dissInt_gas = dissSum_gas*infoStep + dissInt_gas_old;
		dissInt_int = dissSum_int*infoStep + dissInt_int_old;
		
		diffDissInt_liq = diffDissSum_liq*infoStep + diffDissInt_liq_old;
		diffDissInt_gas = diffDissSum_gas*infoStep + diffDissInt_gas_old;
		diffDissInt_int = diffDissSum_int*infoStep + diffDissInt_int_old;
				
		energyDiff = (energy_n - energy_n_old)/infoStep; //difference in total free energy
		kinDiff = (kin - kin_old)/infoStep;
		//kinGDiff = (kinG - kinG_old)/infoStep;
		
		if (t==0){
			kinDiff = 0;
		}
		
		/*if (t == equilTime) {
			energy_init = energy;
			kin_init = kin;
		}*/
		
		/*String fileName3("dissipation.dat");
		ofstream file3(fileName3.get(), ios::app);
		file3.precision(7);
		//file3<< t <<" " << RXcm<<" "<<RVXcm<<" "<<RXcm1+0.4244132*dropletR*sqrt2 <<" "<<RVXcm1<<" "<<RVZcm<<" "<< ECCX << " "<<kin<<" "<<" "<<kin_CM<<endl; 
		//file3<< t << "   " << dissSum_liq << " " << dissSum_int<< " " << dissSum_gas <<"  "<< dissSum_liq/mLiq << " " << dissSum_int/mInt<< " " << dissSum_gas/mGas << "     " << dissInt_liq << " " << dissInt_int<< " " << dissInt_gas << "  " << dissInt_liq/mLiq << " " << dissInt_int/mInt<< " " << dissInt_gas/mGas << "           " << energy_n << " " << kin << " " << kin_gas << " " << energy_n_init << "     " << energyDiff <<" "<< kinDiff << " " << kin_gasDiff << "     "  << dissSum2 << " "<< dissInt2 << endl; 
		file3<< t << "   " << dissSum_liq << " " << dissSum_int<< " " << dissSum_gas <<"  "<< dissSum_liq/mLiq << " " << dissSum_int/mInt<< " " << dissSum_gas/mGas << "     " << dissInt_liq << " " << dissInt_int<< " " << dissInt_gas << "  " << dissInt_liq/mLiq << " " << dissInt_int/mInt<< " " << dissInt_gas/mGas << "           " << energy_n << " " << energy_n_init << " " << kin_liq << " " << kin_int << " " << kin_gas << " " << kin << " " << kin_init << "     " << energyDiff <<" "<< kinDiff << endl; 

		file3.close();*/
		
		/*String fileName4("diffdiss.dat");
		ofstream file4(fileName4.get(), ios::app);
		file4.precision(7);
		file4<< t << "   " << diffDissSum_liq << " " << diffDissSum_int<< " " << diffDissSum_gas <<"  "<< diffDissSum_liq/mLiq << " " << diffDissSum_int/mInt<< " " << diffDissSum_gas/mGas << "     " << diffDissInt_liq << " " << diffDissInt_int<< " " << diffDissInt_gas << "  " << diffDissInt_liq/mLiq << " " << diffDissInt_int/mInt<< " " << diffDissInt_gas/mGas << "      " << LFree << " " << LFree_init << " " << diffDiss_init << "   " << diffDiss_surf << " " << diss_surf << endl; 		
		file4.close();*/
		
		kinDiff = (kin - kin_old)/infoStep;
				
		dissInt_liq_old  = dissInt_liq;
		dissInt_gas_old  = dissInt_gas;
		dissInt_int_old  = dissInt_int;
		diffDissInt_liq_old  = diffDissInt_liq;
		diffDissInt_gas_old  = diffDissInt_gas;
		diffDissInt_int_old  = diffDissInt_int;
		
		energy_n_old = energy_n;  
		kin_old = kin;
		
		
	}
	
}

