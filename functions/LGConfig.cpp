//LGConfig.cpp: initialise the Binary fluid configuration used in initialise()
//Marco, 3 July 2011

#include "wet.h"

void wet::LGConfig(void) 
{	
	//cout << "Process "<< rank << ": LGconfig started" << endl;
	
		
	double teta; 
	double Rs, Hs, Xs, Rs2;
	int kindN=0;
	
	
	{
						
		areaFrac=1.0;
	//	cout << "LG Config check3 " << endl;
		teta_CB = acos(areaFrac * cos(teta1) + areaFrac - 1.0);
	//	cout << "LG Config check4 " << endl;
		
		if (teta_CB_fix != 0) { 
			teta = teta_CB_fix;
			if (rank==ROOT && k==k1) cout << "fixed Theta_CB to " << " " << teta*180.0/M_PI << endl;
		}
		if (teta_CB_fix == 0) {
			teta = teta_CB;
		}
		//if (teta_CB_fix*180/M_PI <= 179.9) dropletCenterZ = 0;
		
		if (dimensions == 2) {
			Rs  = sqrt( M_PI/(teta    - 0.5*sin(2*teta   )) )*dropletR;
			Rs2 = sqrt( M_PI/(teta_CB - 0.5*sin(2*teta_CB)) )*dropletR;
		}
		if (dimensions == 3) {
			Rs  = pow( 4.0/(2.0-3.0*cos(teta)   + pow(cos(teta   ), 3.0) )  ,(1.0/3.0) ) *dropletR;
			Rs2 = pow( 4.0/(2.0-3.0*cos(teta_CB)+ pow(cos(teta_CB), 3.0) )  ,(1.0/3.0) ) *dropletR;
		}
		
		
		Hs = Rs* cos(M_PI-teta) + Dh + dropletCenterZ;  // Dh is height of floor/posts, but from rectposts.dat
		if (teta_CB_fix*180/M_PI <= 179.9)
			Hs = Rs* cos(M_PI-teta) + Dh;
		
		Xs = dropletCenterX - (Rs2 -dropletR) ;
		
		//cout << "Process " << rank <<": LG Config check5 " << endl;
		
				
	}

	computeCoordinate(k);
	
	
		
	
	if ((xk-dropletCenterX)*(xk-dropletCenterX)+(yk-dropletCenterY)*(yk-dropletCenterY)<=dropletR*dropletR) {
        kindN=1;
    }
	

	if(kindN==1 )
	{
		n[k] = 1.0;						//TOTAL density at node k for fluid f
		p[k] = 1.0;						//order parameter 
		uxs[k] = initUX; 
		uys[k] = initUY; 
		uzs[k] = initUZ;	            //Velocity at node k for fluid 1
	} 

	if(kindN==0) 
	{
		n[k] = 1.0;
		p[k] = -1.0;
		uxs[k] = initUX;
		uys[k] = initUY; 
		uzs[k] = initUZ;	  
	} 
	//cout << "Process " << rank <<": LG Config check 7 " << endl;
    //cout << p[k] << endl;
	if (strcmp(geometry,"TESTSYSTEM") == 0) {
		
		n[k] = 1.0;
		p[k] = 1.0;
		
		computeCoordinate(k);
		uxs[k] = initUX*( -1.0* fabs(zk - (LZ-Dh)/2.0)/(LZ -Dh)*2 + 0.5);
		if (zk < Dh) {
			uxs[k] = 0.0;
		}
		//uxs[k] = initUX;
		//cout << k << " " << uxs[k] << endl;
		uys[k] = 0.0;
		uzs[k] = 0.0;
	}
	if (strcmp(geometry,"TESTSYSTEMb") == 0) {
		
		n[k] = 1.0;
		p[k] = 1.0;
		
		computeCoordinate(k);
		uxs[k] = initUX*((zk-Dh)/(LZ-Dh));
		if (zk < Dh) {
			uxs[k] = 0.0;
		}
		uys[k] = 0.0;
		uzs[k] = 0.0;
	}
	if (strcmp(geometry,"TESTSYSTEM2") == 0) {
		
		n[k] = 1.0;
		p[k] = -1.0;
		
		computeCoordinate(k);
		if ( (zk > 1*LZ/8 && zk < 3*LZ/8) || (zk > 5*LZ/8 && zk < 7*LZ/8) ) {
			n[k] = 1.0;
			p[k] = 1.0;
			
		}
		uxs[k] = initUX*( -1.0* fabs(zk - (LZ-Dh)/2.0)/(LZ -Dh)*2 + 0.5);
		if (zk < Dh) {
			uxs[k] = 0.0L;
		}
		
		uys[k] = 0.0L;
		uzs[k] = 0.0L;
		
	}
	if (strcmp(geometry,"TESTSYSTEM2eq") == 0) {
		
		n[k] = 1.0;
		p[k] = -1.0;
		
		computeCoordinate(k);
		if ((zk > 1*LZ/8 && zk < 3*LZ/8) || (zk > 5*LZ/8 && zk < 7*LZ/8)) {
			n[k] = 1.0;
			p[k] = 1.0;
			
		}
		
		uxs[k] = 0.0L;
		uys[k] = 0.0L;
		uzs[k] = 0.0L;
		
	}
	
    if (strcmp(geometry,"TESTSYSTEM3") == 0) {
		
		n[k] = 1.0;
		p[k] = -1.0;
		
		computeCoordinate(k);
		if (xk > LX/3 && xk < 2*LX/3) {
			n[k] = 1.0;
			p[k] = 1.0;
			
		}
		//uxs[k] = initUX*( -1.0* fabs(zk - (LZ-Dh)/2.0)/(LZ -Dh)*2 + 0.5);
		//if (zk < Dh) {
        uxs[k] = 0.0;
		//}
		
		uys[k] = 0.0;
		uzs[k] = 0.0;
	}
	
	if (strcmp(geometry,"TESTSYSTEM3eq") == 0) {
		
		//cout << "Process " << rank <<": LG Config check 8 " << endl;
		
		n[k] = 1.0;
		p[k] = -1.0;
		
		computeCoordinate(k);
		if (xk > 2*LX/8 && xk < 6*LX/8) {
			n[k] = 1.0;
			p[k] = 1.0;
			
		}
		uxs[k] = 0.0;
		uys[k] = 0.0;
		uzs[k] = 0.0;
		
		//cout << "Process " << rank <<": LG Config check 9 " << endl;
		
	}
	if (strcmp(geometry,"TESTSYSTEM4") == 0) {
		
		n[k] = 1.0;
		p[k] = -1.0;
		
		computeCoordinate(k);
		if((xk-dropletCenterX)*(xk-dropletCenterX) + (zk-dropletCenterZ)*(zk-dropletCenterZ) - Rs*Rs <= 0 && zk>=Dh ) //Dx, Dy post parameters, hijacked here for ellipse so i don't have to introduce new parameters
		{
			n[k] = 1.0;
			p[k] = 1.0;
			
			uxs[k] = initUX;
			uys[k] = initUY;
			uzs[k] = initUZ;
		}
		
		
	}
	if (strcmp(geometry,"TESTSYSTEM4eq") == 0) {
		
		n[k] = 1.0;
		p[k] = -1.0;
		uxs[k] = 0.0;
		uys[k] = 0.0;
		uzs[k] = 0.0;
		
		
		computeCoordinate(k);
		if((xk-dropletCenterX)*(xk-dropletCenterX) + (zk-dropletCenterZ)*(zk-dropletCenterZ) - Rs*Rs <= 0 && zk>=Dh )
		{
			n[k] = 1.0;
			p[k] = 1.0;
			
		}
		
		
	}
	if (strcmp(geometry,"TESTSYSTEM5") == 0) {
		
		n[k] = 1.0;
		p[k] = -1.0;
		
		//cout << "*  *  * Process " << rank << ": , t = " << t << " LGconfig TS5 " << endl;
		computeCoordinate(k);
		if((xk-dropletCenterX)*(xk-dropletCenterX)/(Dx*Dx) + (zk-dropletCenterZ)*(zk-dropletCenterZ)/(Dy*Dy) - 1.0 <= 0 && zk>=Dh ) //Dx, Dy post parameters, hijacked here for ellipse so i don't have to introduce new parameters
		{
			n[k] = 1.0;
			p[k] = 1.0;
			
			uxs[k] = initUX;
			uys[k] = initUY;
			uzs[k] = initUZ;
		}
		
		
	}
	if (strcmp(geometry,"TESTSYSTEM5eq") == 0) {
		
		//cout << "*  *  * Process " << rank << ": , t = " << t << " LGconfig TS5eq " << endl;
		
		n[k] = 1.0;
		p[k] = -1.0;
		uxs[k] = 0.0L;
		uys[k] = 0.0L;
		uzs[k] = 0.0L;
		
		
		computeCoordinate(k);
		if((xk-dropletCenterX)*(xk-dropletCenterX) + (zk-dropletCenterZ)*(zk-dropletCenterZ) - Rs*Rs <= 0 && zk>=Dh ) //Dx, Dy post parameters, hijacked here for ellipse so i don't have to introduce new parameters
		{
			n[k] = 1.0;
			p[k] = 1.0;
			
			uxs[k] = 0.0L;
			uys[k] = 0.0L;
			uzs[k] = 0.0L;
		}
		
		
	}
	if (strcmp(geometry,"POISEUILLE") == 0) {
		
		n[k] = 1.0;
		p[k] = 1.0;
		
		computeCoordinate(k);
		//uxs[k] = initUX*( -1.0* fabs(zk - (LZ-Dh)/2.0)/(LZ -Dh)*2 ); //we don't know what external force results in what v_max
		uxs[k] = 0.0;
		uys[k] = 0.0;
		uzs[k] = 0.0;
	}
	//cout << "Process " << rank <<": k = " << k << " LG Config check 10 " << endl;
	
	}
