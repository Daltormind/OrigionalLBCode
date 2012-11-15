//LGConfig.cpp: initialise the Binary fluid configuration used in initialise()
//Marco, 3 July 2011

#include "wet.h"

void wet::LGConfigRev(void) 
{	
	//cout << "Process "<< rank << ": LGconfig started" << endl;
	
	//cout << "LG Config check1 " << endl;
	
	double teta; 
	double Rs, Hs, Xs, Rs2;
	int kindN=0;
	//cout << "LG Config check2 " << endl;
	
	{
						
		if (strcmp(geometry,"FLOOR") == 0) {
			areaFrac = 1.0;
		}
		
		if (strcmp(geometry,"YRIDGES") == 0) {
			areaFrac = double(Dy)/double(PeriodY);
		}
		
		if (strcmp(geometry,"XRIDGES") == 0) {
			areaFrac = double(Dx)/double(PeriodX);
		}
		
		if (strcmp(geometry,"POSTS") == 0) {
			areaFrac = (double(Dx)/double(PeriodX))*(double(Dy)/double(PeriodY)) ;
		}
		
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
		if (teta_CB_fix*180/M_PI <= 179.9) dropletCenterZ = 0;
		
		Rs = sqrt( M_PI/(teta-0.5*sin(2*teta)) )*dropletR;
		Rs2= sqrt( M_PI/(teta_CB-0.5*sin(2*teta_CB)) )*dropletR;	
		Hs = Rs* cos(M_PI-teta) + Dh + dropletCenterZ;  // Dh is height of floor/posts, but from rectposts.dat
		Xs = dropletCenterX - (Rs2 -dropletR) ; 
		
		//cout << "Process " << rank <<": LG Config check5 " << endl;
		
				
	}

	computeCoordinate(k);
	
	
	if (dimensions == 3) {
		
		//Create a drop of kindN=1, resting on the floor
		if((xk-Xs)*(xk-Xs) + (yk-dropletCenterY)*(yk-dropletCenterY)+ (zk-Hs)*(zk-Hs) - Rs*Rs <= 0 && zk>Dh ) 
			kindN =1;
		//Periodical condition on y axis
		if((xk-Xs)*(xk-Xs) + (LY-yk+dropletCenterY)*(LY-yk+dropletCenterY) + (zk-Hs)*(zk-Hs) - Rs*Rs <= 0 && zk>Dh) 
			kindN=1;
		if((xk-Xs)*(xk-Xs) + (yk+LY-dropletCenterY)*(yk+LY-dropletCenterY) + (zk-Hs)*(zk-Hs) - Rs*Rs <= 0 && zk>Dh) 
			kindN=1;
		//Periodical condition on x axis        
		if((LX-xk+Xs)*(LX-xk+Xs) + (yk-dropletCenterY)*(yk-dropletCenterY) + (zk-Hs)*(zk-Hs) - Rs*Rs <= 0 && zk>Dh) 
			kindN=1;
		if((xk+LX-Xs)*(xk+LX-Xs) + (yk-dropletCenterY)*(yk-dropletCenterY) + (zk-Hs)*(zk-Hs) - Rs*Rs <= 0 && zk>Dh) 
			kindN=1;
		//Periodical condition on both x and y axis
		if((LX-xk+Xs)*(LX-xk+Xs)+(LY-yk+dropletCenterY)*(LY-yk+dropletCenterY)+(zk-Hs)*(zk-Hs)-Rs*Rs <= 0 && zk>Dh) 
			kindN=1;
		if((xk+LX-Xs)*(xk+LX-Xs)+(yk+LY-dropletCenterY)*(yk+LY-dropletCenterY)+(zk-Hs)*(zk-Hs)-Rs*Rs <= 0 && zk>Dh) 
			kindN=1;
		if((xk+LX-Xs)*(xk+LX-Xs)+(LY-yk+dropletCenterY)*(LY-yk+dropletCenterY)+(zk-Hs)*(zk-Hs)-Rs*Rs <= 0 && zk>Dh) 
			kindN=1;
		if((LX-xk+Xs)*(LX-xk+Xs)+(yk+LY-dropletCenterY)*(yk+LY-dropletCenterY)+(zk-Hs)*(zk-Hs)-Rs*Rs <= 0 && zk>Dh) 
			kindN=1;

	
		if (equilfirst == false) {  //create second drop 
			
			//Create a drop of kindN=1, resting on the floor
			if((xk-(LX-Xs))*(xk-(LX-Xs)) +(yk-dropletCenterY)*(yk-dropletCenterY)+ (zk-Hs)*(zk-Hs) - Rs*Rs <= 0 && zk>Dh ) 
				kindN =1;
			//Periodical condition on y axis
			if((xk-(LX-Xs))*(xk-(LX-Xs)) + (LY-yk+dropletCenterY)*(LY-yk+dropletCenterY) + (zk-Hs)*(zk-Hs) - Rs*Rs <= 0 && zk>Dh) 
				kindN=1;
			if((xk-(LX-Xs))*(xk-(LX-Xs)) + (yk+LY-dropletCenterY)*(yk+LY-dropletCenterY) + (zk-Hs)*(zk-Hs) - Rs*Rs <= 0 && zk>Dh) 
				kindN=1;
			//Periodical condition on x axis        
			if((LX-xk+(LX-Xs))*(LX-xk+(LX-Xs)) + (yk-dropletCenterY)*(yk-dropletCenterY) + (zk-Hs)*(zk-Hs) - Rs*Rs <= 0 && zk>Dh) 
				kindN=1;
			if((xk+LX-(LX-Xs))*(xk+LX-(LX-Xs)) + (yk-dropletCenterY)*(yk-dropletCenterY) + (zk-Hs)*(zk-Hs) - Rs*Rs <= 0 && zk>Dh) 
				kindN=1;
			//Periodical condition on both x and y axis
			if((LX-xk+(LX-Xs))*(LX-xk+(LX-Xs)) + (LY-yk+dropletCenterY)*(LY-yk+dropletCenterY) + (zk-Hs)*(zk-Hs) - Rs*Rs <= 0 && zk>Dh) 
				kindN=1;
			if((xk+LX-(LX-Xs))*(xk+LX-(LX-Xs)) + (yk+LY-dropletCenterY)*(yk+LY-dropletCenterY) + (zk-Hs)*(zk-Hs) - Rs*Rs <= 0 && zk>Dh) 
				kindN=1;
			if((xk+LX-(LX-Xs))*(xk+LX-(LX-Xs)) + (LY-yk+dropletCenterY)*(LY-yk+dropletCenterY) + (zk-Hs)*(zk-Hs) - Rs*Rs <= 0 && zk>Dh) 
				kindN=1;
			if((LX-xk+(LX-Xs))*(LX-xk+(LX-Xs)) + (yk+LY-dropletCenterY)*(yk+LY-dropletCenterY) + (zk-Hs)*(zk-Hs) - Rs*Rs <= 0 && zk>Dh) 
				kindN=1;
			
		}
	 
	}
	//cout << "Process " << rank <<": LG Config check6 " << endl;	
	
	if (dimensions == 2){
		
		//Create a drop of kindN=1, resting on the floor
		if((xk-Xs)*(xk-Xs) + (zk-Hs)*(zk-Hs) - Rs*Rs <= 0 && zk>Dh ) 
			kindN =1;
		//Periodical condition on y axis
		if((xk-Xs)*(xk-Xs) + (zk-Hs)*(zk-Hs) - Rs*Rs <= 0 && zk>Dh) 
			kindN=1;
		if((xk-Xs)*(xk-Xs) + (zk-Hs)*(zk-Hs) - Rs*Rs <= 0 && zk>Dh) 
			kindN=1;
		//Periodical condition on x axis        
		if((LX-xk+Xs)*(LX-xk+Xs) + (zk-Hs)*(zk-Hs) - Rs*Rs <= 0 && zk>Dh) 
			kindN=1;
		if((xk+LX-Xs)*(xk+LX-Xs) + (zk-Hs)*(zk-Hs) - Rs*Rs <= 0 && zk>Dh) 
			kindN=1;
		//Periodical condition on both x and y axis
		if((LX-xk+Xs)*(LX-xk+Xs) + (zk-Hs)*(zk-Hs) - Rs*Rs <= 0 && zk>Dh) 
			kindN=1;
		if((xk+LX-Xs)*(xk+LX-Xs) + (zk-Hs)*(zk-Hs) - Rs*Rs <= 0 && zk>Dh) 
			kindN=1;
		if((xk+LX-Xs)*(xk+LX-Xs) + (zk-Hs)*(zk-Hs) - Rs*Rs <= 0 && zk>Dh) 
			kindN=1;
		if((LX-xk+Xs)*(LX-xk+Xs) + (zk-Hs)*(zk-Hs) - Rs*Rs <= 0 && zk>Dh) 
			kindN=1;
		
		
		if (equilfirst == false) {  //create second drop 
			
			//Create a drop of kindN=1, resting on the floor
			if((xk-(LX-Xs))*(xk-(LX-Xs)) + (zk-Hs)*(zk-Hs) - Rs*Rs <= 0 && zk>Dh ) 
				kindN =1;
			//Periodical condition on y axis
			if((xk-(LX-Xs))*(xk-(LX-Xs)) + (zk-Hs)*(zk-Hs) - Rs*Rs <= 0 && zk>Dh) 
				kindN=1;
			if((xk-(LX-Xs))*(xk-(LX-Xs)) + (zk-Hs)*(zk-Hs) - Rs*Rs <= 0 && zk>Dh) 
				kindN=1;
			//Periodical condition on x axis        
			if((LX-xk+(LX-Xs))*(LX-xk+(LX-Xs)) + (zk-Hs)*(zk-Hs) - Rs*Rs <= 0 && zk>Dh) 
				kindN=1;
			if((xk+LX-(LX-Xs))*(xk+LX-(LX-Xs)) + (zk-Hs)*(zk-Hs) - Rs*Rs <= 0 && zk>Dh) 
				kindN=1;
			//Periodical condition on both x and y axis
			if((LX-xk+(LX-Xs))*(LX-xk+(LX-Xs)) + (zk-Hs)*(zk-Hs) - Rs*Rs <= 0 && zk>Dh) 
				kindN=1;
			if((xk+LX-(LX-Xs))*(xk+LX-(LX-Xs)) + (zk-Hs)*(zk-Hs) - Rs*Rs <= 0 && zk>Dh) 
				kindN=1;
			if((xk+LX-(LX-Xs))*(xk+LX-(LX-Xs)) + (zk-Hs)*(zk-Hs) - Rs*Rs <= 0 && zk>Dh) 
				kindN=1;
			if((LX-xk+(LX-Xs))*(LX-xk+(LX-Xs)) + (zk-Hs)*(zk-Hs) - Rs*Rs <= 0 && zk>Dh) 
				kindN=1;
			
		}
	 
	}
	//cout << "LG Config check7 " << endl;
	

	if(kindN==0 )
	{
		n[k] = 1.0;						//TOTAL density at node k for fluid f
		p[k] = 1.0;						//order parameter 
		uxs[k] = initUX; 
		uys[k] = initUY; 
		uzs[k] = initUZ;	            //Velocity at node k for fluid 1
	} 

	if(kindN==1) 
	{
		n[k] = 1.0;
		p[k] = -1.0;
		uxs[k] = 0.0; 
		uys[k] = 0.0; 
		uzs[k] = 0.0;	  
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
		uxs[k] = initUX*( -1.0* fabs(zk - (LZ-Dh)/2.0)/(LZ -Dh)*2 + 0.5); 
		if (zk < Dh) {
			uxs[k] = 0.0;
		}
		
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
		if((xk-dropletCenterX)*(xk-dropletCenterX) + (zk-dropletCenterZ)*(zk-dropletCenterZ) - Rs*Rs <= 0 && zk>Dh ) //Dx, Dy post parameters, hijacked here for ellipse so i don't have to introduce new parameters
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
		if((xk-dropletCenterX)*(xk-dropletCenterX) + (zk-dropletCenterZ)*(zk-dropletCenterZ) - Rs*Rs <= 0 && zk>Dh ) 
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
		if((xk-dropletCenterX)*(xk-dropletCenterX)/(Dx*Dx) + (zk-dropletCenterZ)*(zk-dropletCenterZ)/(Dy*Dy) - 1.0 <= 0 && zk>Dh ) //Dx, Dy post parameters, hijacked here for ellipse so i don't have to introduce new parameters
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
		if((xk-dropletCenterX)*(xk-dropletCenterX) + (zk-dropletCenterZ)*(zk-dropletCenterZ) - Rs*Rs <= 0 && zk>Dh ) //Dx, Dy post parameters, hijacked here for ellipse so i don't have to introduce new parameters
		{
			n[k] = 1.0;						
			p[k] = 1.0;
			
			uxs[k] = 0.0L; 
			uys[k] = 0.0L; 
			uzs[k] = 0.0L;	
		}		
		
		
	}
	
	//cout << "Process " << rank <<": k = " << k << " LG Config check 10 " << endl;
}
