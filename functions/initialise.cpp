//ReadInput.cpp: initialise most of computational variables
//Marco, 28 June 2011

#include "wet.h"

void wet::initialise(void)
{  
	long d, u, r, l, q, w;		//Index that works out the neighbours of k {(l,r),(d,u),(q,w)}
	MPI_Status statusLeft, statusRight;

	//periodic boundary condition in processors
	if(rank==0)
		leftProcess=size-1;
	else
		leftProcess=rank-1;
	
	if(rank==size-1)
		rightProcess=0;
	else
		rightProcess=rank+1;
	cout << "Initialise: Process " << rank << " : rightProcess: " << rightProcess << " leftProcess: " << leftProcess << endl;
	
	
		
		LX = LXc;
		teta1 = tetac;
	
	
	
	double alpha=acos(sin(teta1)*sin(teta1));
	phi11 = -2.0*sign(M_PI/2.0L-teta1)*sqrt(A/2/kappa_p)*sqrt(cos(alpha/3)*(1-cos(alpha/3)));  
	cout << "sign " << sign(M_PI/2.0L-teta1) << endl;
    cout << "phi11 init=" << phi11 << endl;
    cout << "alpha=" << alpha << endl;
    alpha=acos(sin(teta2)*sin(teta2));
	phi12 = -2.0*sign(M_PI/2.0L-teta2)*sqrt(A/2/kappa_p)*sqrt(cos(alpha/3)*(1-cos(alpha/3))); 
	
	p_thresh_n = tanh(1.0);
	//MPI_Barrier(MPI_COMM_WORLD);
	
	cout << "Process "<< rank <<": LX" << LX << " , LXc " << LXc << endl; 

	
	
	N = LX*LY*LZ;			//toal number of lattice points in the simulation

	if(size>LX)
	{
		cout << "Process " << rank << ": ERROR, too many processes for this lattice. Abort!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 112);
	}

	if(size > (int)(LX/3))
		cout << "Process "<< rank <<": WARNING, very large number of processes for this lattice" << endl; 

	if(rank==ROOT)
		processN=LY*LZ*((LX-LX%size)/size + LX%size + 2);	//Number of YZ plane for root process plus neighboors
	else
		processN=LY*LZ*((LX-LX%size)/size+2);			    //Number of YZ plane for other processes

	k1 = LY*LZ;			//k where real lattice starts
	k2 = processN-LY*LZ;		//k where real lattice ends

	cout << "Process " << rank << ": k1 = " << k1 << endl;
	cout << "Process " << rank << ": k2 = " << k2 << endl;

	dx = 1.0; dt = 1.0 , mo=1.0;		//initialization of step and time step lattice, used in equilibrium function calculation, in equilibrium()
	c = dx/dt; c2 = c*c;		//velocity and its square, used in equilibrium function calculation, in equilibrium()

	z1=1.0/(6*c);			//Costant used in equilibrium function calculation, in equilibrium()
	z2=1.0/(12*c2);			//Costant used in equilibrium function calculation, in equilibrium()
	z3=5.0/(12*c2);			//Costant used in equilibrium function calculation, in equilibrium()
	z4=1.0/(3*c2);			//Costant used in equilibrium function calculation, in equilibrium()
	z5=1.0/(12*c);			//Costant used in equilibrium function calculation, in equilibrium()
	z6=1.0/(24*c2);			//Costant used in equilibrium function calculation, in equilibrium()
	z7=1.0/(4*c2);			//Costant used in equilibrium function calculation, in equilibrium()

	cout << "Process "<< rank << ": allocation memory for " <<processN <<" lattice node..." <<endl;	
	// Allocate memory for variables                           
	//ff(i) is a pointer to a double vector of dimension N for particle of fluid f, with direction i, after streaming and collision
	ff0 = new double[processN]; 
	ff1 = new double[processN]; ff2 = new double[processN]; ff3 = new double[processN]; 
	ff4 = new double[processN]; ff5 = new double[processN]; ff6 = new double[processN]; 
	ffa = new double[processN]; ffb = new double[processN]; ffc = new double[processN]; 
	ffd = new double[processN]; ffe = new double[processN]; fff = new double[processN];
	ffg = new double[processN]; ffh = new double[processN]; ffi = new double[processN];
	ffj = new double[processN]; ffk = new double[processN]; ffl = new double[processN];

	//gg(i) is a pointer to a double vector of dimension N for particle of fluid g, with direction i, after streaming and collision
	gg0 = new double[processN]; 
	gg1 = new double[processN]; gg2 = new double[processN]; gg3 = new double[processN]; 
	gg4 = new double[processN]; gg5 = new double[processN]; gg6 = new double[processN]; 
	gga = new double[processN]; ggb = new double[processN]; ggc = new double[processN]; 
	ggd = new double[processN]; gge = new double[processN]; ggf = new double[processN];
	ggg = new double[processN]; ggh = new double[processN]; ggi = new double[processN];
	ggj = new double[processN]; ggk = new double[processN]; ggl = new double[processN];

	//fn(i) is a pointer to a double vector of dimension N for particle of fluid f, with direction i, after collision, before streaming
	fn0 = new double[processN]; 
	fn1 = new double[processN]; fn2 = new double[processN]; fn3 = new double[processN]; 
	fn4 = new double[processN]; fn5 = new double[processN]; fn6 = new double[processN]; 
	fna = new double[processN]; fnb = new double[processN]; fnc = new double[processN]; 
	fnd = new double[processN]; fne = new double[processN]; fnf = new double[processN];
	fng = new double[processN]; fnh = new double[processN]; fni = new double[processN];
	fnj = new double[processN]; fnk = new double[processN]; fnl = new double[processN];
	
	//gn(i) is a pointer to a double vector of dimension N for particle of fluid f, with direction i, after collision, before streaming
	gn0 = new double[processN]; 
	gn1 = new double[processN]; gn2 = new double[processN]; gn3 = new double[processN]; 
	gn4 = new double[processN]; gn5 = new double[processN]; gn6 = new double[processN]; 
	gna = new double[processN]; gnb = new double[processN]; gnc = new double[processN]; 
	gnd = new double[processN]; gne = new double[processN]; gnf = new double[processN];
	gng = new double[processN]; gnh = new double[processN]; gni = new double[processN];
	gnj = new double[processN]; gnk = new double[processN]; gnl = new double[processN];

  
	n = new double[processN];			//Total density of lattice points	
	p = new double[processN];
	//Total density of phi
	/*
	if(afterequilflag == false){
		ndup = new double[processN];
		pdup = new double[processN];
	}
	*/
	
	uxs = new double[processN]; 	//Velocity of lattice points			(SISTEMA)
	uys = new double[processN]; 
	uzs = new double[processN]; 
  	
	//Pointers to a vector of dimendion N that represent forces acting on lattice points, in initialise()	      (SISTEMA)
	forcex = new double[processN]; forcey = new double[processN]; forcez = new double[processN];
	//Pointers to a vector of dimension N that is the neighbours of each lattice point, in initialise()
	dd1 = new long[processN]; dd2 = new long[processN]; dd3 = new long[processN]; 
	dd4 = new long[processN]; dd5 = new long[processN]; dd6 = new long[processN]; 
	dda = new long[processN]; ddb = new long[processN]; ddc = new long[processN]; 
	ddd = new long[processN]; dde = new long[processN]; ddf = new long[processN];
	ddg = new long[processN]; ddh = new long[processN]; ddi = new long[processN];
	ddj = new long[processN]; ddk = new long[processN]; ddl = new long[processN];
  	//Pointers to a vector of dimension N that is lattice ID for topolocially patterned substrate, defined in initialise()
	mask = new int[processN];
	
	
		pGlobal = new double[N];
		nGlobal = new double[N];
		maskGlobal = new int[N];
	
	
	uGlobal = new double[N];
	


	cout << "Process "<< rank << ": allocation completed." << endl;

	cout << "Process " << rank << ": gas and surface initialization..." << endl;

	//Initialise this value to used in equilibrium()
	dp_dx = 0.0; dp_dy = 0.0; dp_dz = 0.0; del2p = 0.0;
	
	//Initialise dissipation variables used in dissipation(), this is for dissipation through shear
	dissipation = new double[processN];
	dissGlobal = new double[N];

	//cout << "Process " << rank << ": , t = " << t << " initialise check 0, " << endl;
	
	//dissipation through diffusion at the interface:
	mu = new double[processN];
	muGlobal = new double[N];
	
	//cout << "Process " << rank << ": , t = " << t << " initialise check 1, "  << endl;
	diffDiss = new double[processN];
	//cout << "Process " << rank << ": , t = " << t << " initialise check 2, "  << endl;
	diffDissGlobal = new double[N];
	/*
	eBal = new double[processN];
	eBal_old = new double[processN];
	eBalDiff = new double[processN];
	eBalDiffGlobal = new double[N];
	*/
	//cout << "Process " << rank << ": , t = " << t << " initialise check 3, "  << endl;
	
	dissSum_liq = 0.0;
	dissSum_gas = 0.0;
	dissSum_int = 0.0;

	dissInt_liq = 0.0;
	dissInt_gas = 0.0;
	dissInt_int = 0.0;

	dissInt_liq_old = 0.0;
	dissInt_gas_old = 0.0;
	dissInt_int_old = 0.0;

	
	diffDissSum_liq = 0.0;
	diffDissSum_gas = 0.0;
	diffDissSum_int = 0.0;
	
	diffDissInt_liq = 0.0;
	diffDissInt_gas = 0.0;
	diffDissInt_int = 0.0;

	diffDissInt_liq_old = 0.0;
	diffDissInt_gas_old = 0.0;
	diffDissInt_int_old = 0.0;

	energy_n_old = 0.0;
	kin_old = 0.0;
		

	//cout << "Process " << rank << ": , t = " << t << " initialise check 4, " << "mu[k1] "<< mu[k1] << endl;
	//cout << "Process " << rank << ": , t = " << t << " initialise check 4, " << "G[0] "<< G[0] << endl;
	//Define the mask for each lattice point, configure surface, configure fluids and define the neighbours
	for (k = k1; k < k2; k++) 
	{	
		mask[k] = 0;
		   
		forcex[k] = G[0]; 
		forcey[k] = G[1]; 
		forcez[k] = G[2];
//cout << "Process " << rank << ": check2 " << endl;
		computeCoordinate(k); 	//calculate xk, yk, zk
		
		//Define the neighbours
		int xkProcess= (int) (k/(float)(LY*LZ));
		l = xkProcess - 1; 
		r = xkProcess + 1; 
		
		if (yk == 0) 
			d = LY - 1;
		else 
			d = yk - 1; 
		if (yk == LY -1) 
			u = 0;
		else 
			u = yk + 1; 
		
		if (zk == 0) 
			q = LZ - 1;
		else 
			q = zk - 1; 
		if (zk == LZ - 1) 
			w = 0;
		else 
			w = zk + 1;        
		
		dd1[k] = zk + yk*LZ + r*LZ*LY; 	  
		dd2[k] = zk + yk*LZ + l*LZ*LY; 
		dd3[k] = zk + u*LZ + xkProcess*LZ*LY; 
		dd4[k] = zk + d*LZ + xkProcess*LZ*LY; 
		dd5[k] = w + yk*LZ + xkProcess*LZ*LY; 
		dd6[k] = q + yk*LZ + xkProcess*LZ*LY; 
		
		dda[k] = zk + u*LZ + r*LZ*LY; 
		ddb[k] = zk + u*LZ + l*LZ*LY; 
		
		ddc[k] = zk + d*LZ + r*LZ*LY; 
		ddd[k] = zk + d*LZ + l*LZ*LY; 
		dde[k] = w + u*LZ + xkProcess*LZ*LY;
		ddf[k] = w + d*LZ + xkProcess*LZ*LY;
		ddg[k] = q + u*LZ + xkProcess*LZ*LY;
		ddh[k] = q + d*LZ + xkProcess*LZ*LY;
		ddi[k] = w + yk*LZ + r*LZ*LY;
		ddj[k] = w + yk*LZ + l*LZ*LY;
		ddk[k] = q + yk*LZ + r*LZ*LY;
		ddl[k] = q + yk*LZ + l*LZ*LY;
		
		
	}
	
	
			
//cout << "Process " << rank << ": check4 " << endl;
		
		
			
			for (k = k1; k < k2; k++){
				
				LGConfig();		//initialise density of two fluids, using the position xk, yk, zk
				//LGConfigRev();
				//if (k==k1) cout << "LG config done" << endl;
				initialiseSurface(); //initialise surface drawing mask=28 if the point xk, yk, zk is part of solid surface
				//cout << "mask" << k << " = " << mask[k];
				if (mask[k] == 28)
				{
					n[k] = 1.0; p[k] = 0.0,uxs[k]=0.0,uys[k]=0.0,uzs[k]=0.0;
				}
			}
			
			exchangeMask();
			
			for (k = k1; k < k2; k++)
				relabel();		//Redrawing mask to put value 1 near the surface
			
			exchangeMask();
			
			
			for (k = k1; k < k2; k++) 
				leakSearch();		//Search for missing point in the relabel function 
			
			generateGlobalMask();
			
		
		/*if( rank==ROOT)
        {for(int i = 0 ; i < LX ; i++) {
            for(int j = 0 ; j < LY ; j++)
            {
                k = j*LZ + i*LY*LZ;
                cout << maskGlobal[k];
                
            }
            cout << endl;
        }}*/
			

			
	cout << "Process " << rank << ": gas and surface initialised." << endl;


	//cout << "------  Process " << rank << ": exchangePhi in initialise..." << endl;
	exchangePhi();	
	//cout << "------  Process " << rank << ": exchangePhi in initialise done" << endl;
	
	//cout << "Process " << rank << ": distribution function initialization..." << endl;

	//Initialise the partial distribution function
	if(afterequilflag == false){ cout << "Process " << rank << ": entered equiliberium Loop" << endl;
		for(k = k1; k < k2; k++)
		{
			/*if (mask[k] == 28) 
			{
				n[k] = 1.0; p[k] = 0.0;
			}
			*/
			
			//cout << "* * * Process " << rank << ", k = " << k << " pre equil() in initialise function initialization..." << endl;
			equilibrium();	//Create equilibrium function
			//cout << "* * * Process " << rank << ", k = " << k << " equil() in initialise function initialization done" << endl;
			//initialise density function of fluid f with equilibrium function
			ff0[k] = fe0; 
			ff1[k] = fe1; ff2[k] = fe2; ff3[k] = fe3; ff4[k] = fe4; ff5[k] = fe5; ff6[k] = fe6; 
			ffa[k] = fea; ffb[k] = feb; ffc[k] = fec; ffd[k] = fed; ffe[k] = fee; fff[k] = fef;
			ffg[k] = feg; ffh[k] = feh; ffi[k] = fei; ffj[k] = fej; ffk[k] = fek; ffl[k] = fel;
			//initialise density function of fluid g with equilibrium function
			gg0[k] = ge0;
			gg1[k] = ge1; gg2[k] = ge2; gg3[k] = ge3; gg4[k] = ge4; gg5[k] = ge5; gg6[k] = ge6; 
			gga[k] = gea; ggb[k] = geb; ggc[k] = gec; ggd[k] = ged; gge[k] = gee; ggf[k] = gef;
			ggg[k] = geg; ggh[k] = geh; ggi[k] = gei; ggj[k] = gej; ggk[k] = gek; ggl[k] = gel;
			
			//cout << "Process " << rank << ": , t = " << t << " * * * in initialise, ff0[k1+Dh] "<< ff0[k1+Dh] << " , ff0[k1+Dh+LY*LZ] "<< ff0[k1+Dh+LY*LZ] << " , ff0[k1+Dh+2*LY*LZ] "<< ff0[k1+Dh+2*LY*LZ]<< endl;		

			//cout << "* * * Process " << rank << ", k = " << k << " PDF in initialise function initialization done" << endl;
			
		}
	}
	exchangeDensities_ffgg();
	cout << "Process " << rank << ": distribution function initialised." << endl;
	
	makematrix();
	
    Dimensionlessnum();
    writeInfoFile();

	/*if(afterequilflag == true){
		for(k = k1; k < k2; k++)
		{
			computeCoordinate(k);
			cout << "AT END OF INITIALISE *** process " << rank << ": uxs at k=" << k << " is " << uxs[k] << endl;
		}
	}*/
}
