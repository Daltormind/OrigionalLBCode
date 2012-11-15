// duplicateDrop.cpp: does what the name says, it copies the pre-equilibration system of size LXc/2*LY*LZ = Nc/2 
// onto the new, double-sized arrays of size N = Nc = LXc*LY*LZ
// Lisa, 31 July 2011

#include "wet.h"



void wet::duplicateDrop1(void)
{
	
	int kdup ;
	int xkdup, ykdup, zkdup;
	//double *nGlobal;
	
	//nGlobal = new double[N];
	
	generatePhiGlobal();
	generateNGlobal();
	generateGlobalMask();
	
	MPI_Bcast(pGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	MPI_Bcast(nGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	MPI_Bcast(maskGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	
	/*for (k=0; k<N;k++){
		cout << " - Process " << rank <<" k " << k << " pG[k] " << pGlobal[k] << endl;
	}*/
	initialise(); //pGlobal, nGlobal, maskGlobal not initialised 
	
	/*for (k=0; k < N/2 ;k++){
		cout << " - - Process " << rank <<" k " << k << " pG[k] " << pGlobal[k] << endl;
	}*/
	
	cout << "Proces " << rank <<": re-initialised" << endl;
	cout << "Proces " << rank << ": new LX = " << LX << endl;
	
	//saveFiles(); //take OUT
	
	if (rank < size/2){          //first half of the box, drop just copied over
		
		
		for (k=k1 ; k< k2; k++) {
			
			
			kdup = k-k1 + rank * (k2-k1);
			
			xkdup =int(kdup/(float) (LZ*LY));
			ykdup =int((kdup-xkdup*LZ*LY)/(float) LZ);
			zkdup =kdup-xkdup*LZ*LY-ykdup*LZ;
			
			
			p[k] = pGlobal[kdup];
			n[k] = nGlobal[kdup];
			mask[k] = maskGlobal[kdup];
			
			//computeCoordinate(k);
			//cout << "-- Process " << rank <<" k= "<< k << " xk yk zk " << xk << " " << yk << " " << zk << " " << "p[k] = " << p[k] << "; kdup= "<< kdup << " xk yk zk " << xkdup << " " << ykdup << " " << zkdup <<" " << "pG[kdup] = " << pGlobal[kdup] << endl;
			//cout << "--- Process " << rank <<": k "<< k << ", kdup in pG "<< kdup <<endl;
			
		}
		
	}
	if (rank >= size/2){
		
		
		for (int kk=k1 ; kk< k2; kk++) {
			
			kdup = kk;
			computeCoordinate(kk);      
			xk = LX-1 - xk;
			invComputeCoordinate();   // sets k, here value on Global arrays; invCompCooordinate not for parallel setup, so convenient
			
			
			xkdup =int(kdup/(float) (LZ*LY));
			ykdup =int((kdup-xkdup*LZ*LY)/(float) LZ);
			zkdup =kdup-xkdup*LZ*LY-ykdup*LZ;
			xkdup =xkdup+rank*LX/size-1%LX;

			
			p[kdup] = pGlobal[k];
			n[kdup] = nGlobal[k];
			mask[kdup] = maskGlobal[k];
			
			//computeCoordinate(k);
			//cout << "<< Process " << rank <<" k= "<< kdup << " xk yk zk " << xkdup << " " << ykdup << " " << zkdup <<" " << "p[k] = " << p[kdup] <<"; kdup= "<< k << " xk yk zk " << xk << " " << yk << " " << zk << " " << "pG[kdup] = " << pGlobal[k] << endl;
			
			
		}
		
	}
	
	if(equilfirst == true && (t - equilTime)<(1)){
		
		if (afterequilflag == true && strcmp(geometry,"TESTSYSTEM2eq") == 0) {
			//cout << "*** process " << rank << ": after equil, TESTSYSTEM2eq " << endl;
			
			for (k = k1; k < k2; k++){
				computeCoordinate(k);
				uxs[k] = initUX*( -1.0* fabs(zk - (LZ-Dh)/2.0)/(LZ -Dh)*2 + 0.5); 
				//cout << "*** process " << rank << ": uxs at k=" << k << " is " << uxs[k] << endl;
				if (zk < Dh) {
					uxs[k] = 0.0;
				}
				uys[k] = 0.0; 
				uzs[k] = 0.0;
			}
		}
		
		if (afterequilflag == true && strcmp(geometry,"TESTSYSTEM4eq") == 0) {
			
			for (k = k1; k < k2; k++){
				if (p[k] >= 0.0L) {
					
					uxs[k] = initUX;
					uys[k] = initUY; 
					uzs[k] = initUZ;
				}
				else{
					
					uxs[k] = 0.0;
					uys[k] = 0.0; 
					uzs[k] = 0.0;
					
				}
				
			}
		}
		
		if (afterequilflag == true && strcmp(geometry,"TESTSYSTEM3eq") == 0) {
			
			for (k = k1; k < k2; k++){
				computeCoordinate(k);
				uxs[k] = initUX*( -1.0* fabs(zk - (LZ-Dh)/2.0)/(LZ -Dh)*2 + 0.5); 
				if (zk < Dh) {
					uxs[k] = 0.0;
				}
				uys[k] = 0.0; 
				uzs[k] = 0.0;
			}		
		}
		if (afterequilflag == true && strcmp(geometry,"TESTSYSTEM5eq") == 0) {
			
			for (k = k1; k < k2; k++){
				if (p[k] < 0.0L) {       //squeeze drop to form ellipse, does only make sense with duplicationtype=2
					
					computeCoordinate(k);
					/*
					 uxs[k] = 0.0L;
					 uys[k] = 0.0L; 
					 uzs[k] = 0.0L;
					 */
					
					if( (xk-dropletCenterX)*(xk-dropletCenterX) + (zk-dropletCenterZ)*(zk-dropletCenterZ) - (dropletR+5)*(dropletR+5) <= 0 || (xk-(LX-dropletCenterX))*(xk-(LX-dropletCenterX)) + (zk-dropletCenterZ)*(zk-dropletCenterZ) - (dropletR+5)*(dropletR+5)  <= 0 ){
						if ( (zk-dropletCenterZ) > dropletR/sqrt(2) && ( fabs(xk-dropletCenterX) <= dropletR/sqrt(2) || fabs(xk-(LX-dropletCenterX)) <= dropletR/sqrt(2)) ) {
							//uxs[k] = 0.0L;
							//uys[k] = 0.0L; 
							uzs[k] = initUZ;
						}
						if ( -(zk-dropletCenterZ) > dropletR/sqrt(2) && ( fabs(xk-dropletCenterX) <= dropletR/sqrt(2) || fabs(xk-(LX-dropletCenterX)) <= dropletR/sqrt(2))) {
							//uxs[k] = 0.0L;
							//uys[k] = 0.0L; 
							uzs[k] = -initUZ;
						}
						if ( fabs(zk-dropletCenterZ) < dropletR/sqrt(2) ) {
							if ( (xk < dropletCenterX) || (LX/2.0 < xk && xk < (LX-dropletCenterX)) ) {
								uxs[k] = initUX;
								//uys[k] = 0.0L; 
								//uzs[k] = 0.0L;
							}
							if ( (dropletCenterX < xk && xk < LX/2.0) || (xk > (LX-dropletCenterX) ) ) {
								uxs[k] = -initUX;
								//uys[k] = 0.0L; 
								//uzs[k] = 0.0L;
							}
							
						}
					}
					
					
				}
			}
		}
		
		
	}
	
	
		
	// have to set equil f and g again after assiging n[k] and p[k], otherwise it's zero. 
	for(k = k1; k < k2; k++){
		
		if (mask[k] == 28) {n[k] = 1.0; p[k] = 0.0;}
		//nn = n[k]; pp = p[k]; ux = uxs[k]; uy = uys[k]; uz = uzs[k]; 
		equilibrium();
		ff0[k] = fe0; 
		ff1[k] = fe1; ff2[k] = fe2; ff3[k] = fe3; ff4[k] = fe4; ff5[k] = fe5; ff6[k] = fe6; 
		ffa[k] = fea; ffb[k] = feb; ffc[k] = fec; ffd[k] = fed; ffe[k] = fee; fff[k] = fef;
		ffg[k] = feg; ffh[k] = feh; ffi[k] = fei; ffj[k] = fej; ffk[k] = fek; ffl[k] = fel;
		
		gg0[k] = ge0;
		gg1[k] = ge1; gg2[k] = ge2; gg3[k] = ge3; gg4[k] = ge4; gg5[k] = ge5; gg6[k] = ge6; 
		gga[k] = gea; ggb[k] = geb; ggc[k] = gec; ggd[k] = ged; gge[k] = gee; ggf[k] = gef;
		ggg[k] = geg; ggh[k] = geh; ggi[k] = gei; ggj[k] = gej; ggk[k] = gek; ggl[k] = gel;
		
	}
	
	/*
	delete []maskGlobal;
	delete []pGlobal;
	delete []nGlobal;
	*/
	pGlobal = new double[N];	
	maskGlobal = new int[N];
	
	generatePhiGlobal();
	generateGlobalMask();
	
	//saveFiles();  // just for testing -  REMOVE! and instead put back in in LGAlg
	if(rank==ROOT)
	{
		int i, j, h;	
		
		char filename[25];
		sprintf(filename, "./densitydata/dd_dupl%ld.m", t);			//Create a name for file that contain data 
		//sprintf(filename, "./dd%ld.m", t);
		
		ofstream file(filename);
		file.precision(4);
		for( h = 0 ; h < LZ ; h++) 
		{   
			file << "ddup" << t << "(:,:," << h+1 << ")=[" << endl;
			for( i = 0 ; i < LX ; i++) 
			{
				for( j = 0 ; j < LY ; j++) 
				{
					k = h + j*LZ + i*LY*LZ;
					if( maskGlobal[k]==28) 
						file << -2.0 << " ";  
					else
						file << pGlobal[k] << " ";
					
				}
				file << endl;
			}
			file <<"];" << endl;
		}
		file.close();
	}
	
	
}
	
