//equilibrium.cpp: calculates the equilibrium dist function, in initialise(). 
//The equilibrium function following Kusumaatmaja, Yeomans, Lattice Boltzmann Simulations of Wetting and Drop Dynamics,
//Springer Berlin / Heidelberg, http://dx.doi.org/10.1007/978-3-642-12203-3_11
//Marco, 29 June 2011

#include "wet.h"

void wet::equilibrium()
{
	//if (rank == ROOT)  cout << "* * * Process " << rank << ", k = " << k << " equilibrium check 0" << endl;

		  
	if (mask[k] != 28) 
	{			
		dp_dx = ((p[dd1[k]] - p[dd2[k]])/6 + (p[dda[k]]+p[ddc[k]]+p[ddi[k]]+p[ddk[k]]-p[ddb[k]]-p[ddd[k]]-p[ddj[k]]-p[ddl[k]])/12);
		dp_dy = ((p[dd3[k]] - p[dd4[k]])/6 + (p[dda[k]]+p[ddb[k]]+p[dde[k]]+p[ddg[k]]-p[ddc[k]]-p[ddd[k]]-p[ddf[k]]-p[ddh[k]])/12);
		dp_dz = ((p[dd5[k]] - p[dd6[k]])/6 + (p[dde[k]]+p[ddf[k]]+p[ddi[k]]+p[ddj[k]]-p[ddg[k]]-p[ddh[k]]-p[ddk[k]]-p[ddl[k]])/12);
		
		
		del2p = (p[dd1[k]] + p[dd2[k]] + p[dd3[k]] + p[dd4[k]] + p[dd5[k]] + p[dd6[k]] - 6*p[k])/3
			+ (p[dda[k]] + p[ddb[k]] + p[ddc[k]] + p[ddd[k]] + p[dde[k]] + p[ddf[k]] + p[ddg[k]] + p[ddh[k]] + p[ddi[k]] + p[ddj[k]] + p[ddk[k]] + p[ddl[k]] - 12*p[k])/6;
	
		// IS THIS CORRECT? MARCO SAYS Factor *2 in second term, because we double count the diagonal bits; CORRECT, cf. OReilly, Beck 2006 (stencils)
		
		//if (rank == ROOT)  cout << "* * * Process " << rank << ", k = " << k << " equilibrium check 1" << endl;
		
		
		mu[k] = - A*p[k] + A*p[k]*p[k]*p[k] - kappa_p*del2p;
		
		if(t % (infoStep) == 0) {
			//if (rank == ROOT)  cout << "* * * Process " << rank << ", k = " << k << " equilibrium check 1b" << endl;
			computeFreeEnergy();  //have to do it here so that we have the gradient of p for the interface energy calculation
			//if (rank == ROOT)  cout << "* * * Process " << rank << ", k = " << k << " equilibrium check 1c" << endl;
		}
	//There is two possible algorithm to use to calculate equilibrium function:
	//This algorithm is slow, but it is "precise", i.e. it does not sum a piece that after it will subtract
		ux2 = uxs[k]*uxs[k]; 
		uy2 = uys[k]*uys[k]; 
		uz2 = uzs[k]*uzs[k];
//if (rank == ROOT)  cout << "* * * Process " << rank << ", k = " << k << " equilibrium check 1d" << endl;
		uxuy2=ux2+uy2;
		uyuz2=uy2+uz2;
		uxuz2=ux2+uz2;
		
		fbulk1 = (n[k]/3-(A/2)*p[k]*p[k] + (3.0*A/4)*p[k]*p[k]*p[k]*p[k]-kappa_p*p[k]*del2p)/(6*c2);
		fbulk2 = fbulk1/2;
		
		//if (rank == ROOT)  cout << "* * * Process " << rank << ", k = " << k << " equilibrium check 1e" << endl;
		nz1=z1*n[k];
		nz2=z2*n[k];
		nz5=z5*n[k];
		nz6ux2=z6*n[k]*ux2;
		nz6uy2=z6*n[k]*uy2;
		nz6uz2=z6*n[k]*uz2;
		nz7uxuy=z7*n[k]*uxs[k]*uys[k];
		nz7uyuz=z7*n[k]*uys[k]*uzs[k];
		nz7uxuz=z7*n[k]*uxs[k]*uzs[k];
		
		//if (rank == ROOT)  cout << "* * * Process " << rank << ", k = " << k << " equilibrium check 2" << endl;
		
		a1 = kappa_p*dp_dx*dp_dx;
		a2 = kappa_p*dp_dy*dp_dy;
		a3 = kappa_p*dp_dz*dp_dz;
		a4 = kappa_p*( dp_dx*dp_dy);
		a5 = kappa_p*( dp_dy*dp_dz);
		a6 = kappa_p*( dp_dz*dp_dx);	

		a12=a1+a2;
		a23=a2+a3;
		a13=a1+a3;
	
	
		fe1 = fbulk1 +nz1*(uxs[k]+ux2/c) -nz2*uyuz2 + z3*a1 -z4*a23;				 
		fe2 = fbulk1 +nz1*(ux2/c-uxs[k]) -nz2*uyuz2 + z3*a1 -z4*a23;  				
		fe3 = fbulk1 +nz1*(uys[k]+uy2/c) -nz2*uxuz2 + z3*a2 -z4*a13;				
		fe4 = fbulk1 +nz1*(uy2/c-uys[k]) -nz2*uxuz2 + z3*a2 -z4*a13;				
		fe5 = fbulk1 +nz1*(uzs[k]+uz2/c) -nz2*uxuy2 + z3*a3 -z4*a12;				
		fe6 = fbulk1 +nz1*(uz2/c-uzs[k]) -nz2*uxuy2 + z3*a3 -z4*a12;				
		fea = fbulk2 +nz5*(uxs[k]+uys[k]+uxuy2/c) +nz7uxuy -nz6uz2 +z2*a3 -z6*a12 +z7*a4;		
		feb = fbulk2 +nz5*(uys[k]+uxuy2/c-uxs[k]) -nz7uxuy -nz6uz2 +z2*a3 -z6*a12 -z7*a4;	
		fec = fbulk2 +nz5*(uxs[k]+uxuy2/c-uys[k]) -nz7uxuy -nz6uz2 +z2*a3 -z6*a12 -z7*a4;
		fed = fbulk2 +nz5*(uxuy2/c-uxs[k]-uys[k]) +nz7uxuy -nz6uz2 +z2*a3 -z6*a12 +z7*a4;
		fee = fbulk2 +nz5*(uys[k]+uzs[k]+uyuz2/c) +nz7uyuz -nz6ux2 +z2*a1 -z6*a23 +z7*a5;
		fef = fbulk2 +nz5*(uzs[k]+uyuz2/c-uys[k]) -nz7uyuz -nz6ux2 +z2*a1 -z6*a23 -z7*a5;
		feg = fbulk2 +nz5*(uys[k]+uyuz2/c-uzs[k]) -nz7uyuz -nz6ux2 +z2*a1 -z6*a23 -z7*a5;
		feh = fbulk2 +nz5*(uyuz2/c-uys[k]-uzs[k]) +nz7uyuz -nz6ux2 +z2*a1 -z6*a23 +z7*a5;
		fei = fbulk2 +nz5*(uxs[k]+uzs[k]+uxuz2/c) +nz7uxuz -nz6uy2 +z2*a2 -z6*a13 +z7*a6; 
		fej = fbulk2 +nz5*(uzs[k]+uxuz2/c-uxs[k]) -nz7uxuz -nz6uy2 +z2*a2 -z6*a13 -z7*a6;
		fek = fbulk2 +nz5*(uxs[k]+uxuz2/c-uzs[k]) -nz7uxuz -nz6uy2 +z2*a2 -z6*a13 -z7*a6;
		fel = fbulk2 +nz5*(uxuz2/c-uxs[k]-uzs[k]) +nz7uxuz -nz6uy2 +z2*a2 -z6*a13 +z7*a6;
	 
		fe0 = n[k] - fe1 - fe2 - fe3 - fe4 - fe5 - fe6 - fea - feb - fec - fed - fee - fef - feg - feh - fei - 	fej - fek - fel;
	  
		//if (rank == ROOT)  cout << "* * * Process " << rank << ", k = " << k << " equilibrium check 3" << endl;
		
		
		gbulk1=gama*(-A*p[k] + A*p[k]*p[k]*p[k] - kappa_p*del2p)/(6*c2);
		gbulk2=gbulk1/2;
	
		pz1=z1*p[k];
		pz2=z2*p[k];
		pz5=z5*p[k];
		pz6ux2=z6*p[k]*ux2;
		pz6uy2=z6*p[k]*uy2;
		pz6uz2=z6*p[k]*uz2;
		pz7uxuy=z7*p[k]*uxs[k]*uys[k];
		pz7uyuz=z7*p[k]*uys[k]*uzs[k];
		pz7uxuz=z7*p[k]*uxs[k]*uzs[k];		
		
		//if (rank == ROOT) cout  << "* * * Process " << rank << ", k = " << k << " equilibrium check 4" << endl;
		
		ge1 = gbulk1 +pz1*(uxs[k]+ux2/c) -pz2*uyuz2; 
		ge2 = gbulk1 +pz1*(ux2/c-uxs[k]) -pz2*uyuz2;
		ge3 = gbulk1 +pz1*(uys[k]+uy2/c) -pz2*uxuz2;
		ge4 = gbulk1 +pz1*(uy2/c-uys[k]) -pz2*uxuz2;
		ge5 = gbulk1 +pz1*(uzs[k]+uz2/c) -pz2*uxuy2;
		ge6 = gbulk1 +pz1*(uz2/c-uzs[k]) -pz2*uxuy2;	
		gea = gbulk2 +pz5*(uxs[k]+uys[k]+uxuy2/c) +pz7uxuy -pz6uz2;	 
		geb = gbulk2 +pz5*(uys[k]+uxuy2/c-uxs[k]) -pz7uxuy -pz6uz2;	
		gec = gbulk2 +pz5*(uxs[k]+uxuy2/c-uys[k]) -pz7uxuy -pz6uz2;	
		ged = gbulk2 +pz5*(uxuy2/c-uxs[k]-uys[k]) +pz7uxuy -pz6uz2;	
		gee = gbulk2 +pz5*(uys[k]+uzs[k]+uyuz2/c) +pz7uyuz -pz6ux2;	
		gef = gbulk2 +pz5*(uzs[k]+uyuz2/c-uys[k]) -pz7uyuz -pz6ux2;	
		geg = gbulk2 +pz5*(uys[k]+uyuz2/c-uzs[k]) -pz7uyuz -pz6ux2;	
		geh = gbulk2 +pz5*(uyuz2/c-uys[k]-uzs[k]) +pz7uyuz -pz6ux2;	
		gei = gbulk2 +pz5*(uxs[k]+uzs[k]+uxuz2/c) +pz7uxuz -pz6uy2;	 
		gej = gbulk2 +pz5*(uzs[k]+uxuz2/c-uxs[k]) -pz7uxuz -pz6uy2;	
		gek = gbulk2 +pz5*(uxs[k]+uxuz2/c-uzs[k]) -pz7uxuz -pz6uy2;	
		gel = gbulk2 +pz5*(uxuz2/c-uxs[k]-uzs[k]) +pz7uxuz -pz6uy2;	
	  
		ge0 = p[k] - ge1 - ge2 - ge3 - ge4 - ge5 - ge6 - gea - geb - gec - ged - gee - gef - geg - geh - gei - 	gej - gek - gel;
	
		//if (rank == ROOT)  cout << "* * * Process " << rank << ", k = " << k << " equilibrium check 5" << endl;
		
		tau2=(tauliquid+taugas)+p[k]*(tauliquid-taugas);
	
		fx = z2*n[k]*forcex[k]*(1.0-1/tau2);
		fy = z2*n[k]*forcey[k]*(1.0-1/tau2);  
		fz = z2*n[k]*forcez[k]*(1.0-1/tau2);
	     
		fsq = fx*uxs[k]+fy*uys[k]+fz*uzs[k];
		Mxxyy = 3*(fx*uxs[k] + fy*uys[k]);
		Myyzz = 3*(fy*uys[k] + fz*uzs[k]);
		Mxxzz = 3*(fx*uxs[k] + fz*uzs[k]);
		Mxyyx = 3*(fx*uys[k] + fy*uxs[k]);
		Myzzy = 3*(fy*uzs[k] + fz*uys[k]);
		Mxzzx = 3*(fz*uxs[k] + fx*uzs[k]);
	     
		//if (rank == ROOT)  cout << "* * * Process " << rank << ", k = " << k << " equilibrium check 6" << endl;
		
		force1 = 2*(c*fx + 3*fx*uxs[k] - fsq);
		force2 = 2*(3*fx*uxs[k] - fsq - c*fx);
		force3 = 2*(c*fy + 3*fy*uys[k] - fsq);
		force4 = 2*(3*fy*uys[k] - fsq - c*fy);
		force5 = 2*(c*fz + 3*fz*uzs[k] - fsq);
		force6 = 2*(3*fz*uzs[k] - fsq - c*fz);
		forcea = c*(fx + fy) + Mxxyy + Mxyyx - fsq;
		forceb = c*(fy - fx) + Mxxyy - Mxyyx - fsq;
		forcec = c*(fx - fy) + Mxxyy - Mxyyx - fsq;
		forced = c*(-fx -fy) + Mxxyy + Mxyyx - fsq;
		forcee = c*(fy + fz) + Myyzz + Myzzy - fsq;
		forcef = c*(fz - fy) + Myyzz - Myzzy - fsq;
		forceg = c*(fy - fz) + Myyzz - Myzzy - fsq;
		forceh = c*(-fy -fz) + Myyzz + Myzzy - fsq;
		forcei = c*(fz + fx) + Mxxzz + Mxzzx - fsq;
		forcej = c*(fz - fx) + Mxxzz - Mxzzx - fsq;
		forcek = c*(fx - fz) + Mxxzz - Mxzzx - fsq;
		forcel = c*(-fz -fx) + Mxxzz + Mxzzx - fsq;
	
		force0 = -force1-force2-force3-force4-force5-force6-forcea-forceb-forcec-forced-forcee-forcef-forceg-	forceh-forcei-forcej-forcek-forcel ;
	
		//if (rank == ROOT) cout << "* * * Process " << rank << ", k = " << k << " equilibrium check 7" << endl;
		
		//This algorithm is faster, but there is a correlated propagation of errors
	}
	//if (rank == ROOT) cout << "* * * Process " << rank << ", k = " << k << " equilibrium done" << endl;
    }
