//duplicateDrop.cpp: does what the name says (hopefully), by copying all ff, gg, and mask into global arrays first 
// and then duplicating them into double sized arrays
//Lisa, 31 July 2011

#include "wet.h"



void wet::duplicateDrop2()
{
	//copyFFGG(); //generates global arrays, copies ff0..ffl, gg0..ggl into Global arrays
	
	double *ff0Global, *ff1Global, *ff2Global, *ff3Global, *ff4Global, *ff5Global, *ff6Global, *ffaGlobal, *ffbGlobal, *ffcGlobal, *ffdGlobal, *ffeGlobal, *fffGlobal, *ffgGlobal, *ffhGlobal, *ffiGlobal, *ffjGlobal, *ffkGlobal, *fflGlobal;
	double *gg0Global, *gg1Global, *gg2Global, *gg3Global, *gg4Global, *gg5Global, *gg6Global, *ggaGlobal, *ggbGlobal, *ggcGlobal, *ggdGlobal, *ggeGlobal, *ggfGlobal, *gggGlobal, *gghGlobal, *ggiGlobal, *ggjGlobal, *ggkGlobal, *gglGlobal;
	
	//Create global arrays for all the distribution functions
	/*
	
	ff0Global = new double[N];	
	ff0Global = generateGlobal(ff0);
	MPI_Bcast(ff0Global,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	ff0 = duplicateArray(ff0Global);
	delete []ff0Global;
	
	ff1Global = new double[N];	
	ff1Global = generateGlobal(ff1);
	MPI_Bcast(ff1Global,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	ff1 = duplicateArray(ff1Global);
	delete []ff1Global;
	
	ff2Global = new double[N];	
	ff2Global = generateGlobal(ff2);
	MPI_Bcast(ff2Global,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	ff2 = duplicateArray(ff2Global);
	delete []ff2Global;
	
	ff3Global = new double[N];	
	ff3Global = generateGlobal(ff3);
	MPI_Bcast(ff3Global,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	ff3 = duplicateArray(ff3Global);
	delete []ff3Global;
	
	ff4Global = new double[N];	
	ff4Global = generateGlobal(ff4);
	MPI_Bcast(ff4Global,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	ff4 = duplicateArray(ff4Global);
	delete []ff4Global;
	
	ff5Global = new double[N];	
	ff5Global = generateGlobal(ff5);
	MPI_Bcast(ff5Global,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	ff5 = duplicateArray(ff5Global);
	delete []ff5Global;
	
	ff6Global = new double[N];	
	ff6Global = generateGlobal(ff6);
	MPI_Bcast(ff6Global,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	ff6 = duplicateArray(ff6Global);
	delete []ff6Global;
	
	ffaGlobal = new double[N];	
	ffaGlobal = generateGlobal(ffa);
	MPI_Bcast(ffaGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	ffa = duplicateArray(ffaGlobal);
	delete []ffaGlobal;
	
	ffbGlobal = new double[N];	
	ffbGlobal = generateGlobal(ffb);
	MPI_Bcast(ffbGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	ffb = duplicateArray(ffbGlobal);
	delete []ffbGlobal;
	
	ffcGlobal = new double[N];	
	ffcGlobal = generateGlobal(ffc);
	MPI_Bcast(ffcGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	ffc = duplicateArray(ffcGlobal);
	delete []ffcGlobal;
	
	ffdGlobal = new double[N];	
	ffdGlobal = generateGlobal(ffd);
	MPI_Bcast(ffdGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	ffd = duplicateArray(ffdGlobal);
	delete []ffdGlobal;
	
	ffeGlobal = new double[N];	
	ffeGlobal = generateGlobal(ffe);
	MPI_Bcast(ffeGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	ffe = duplicateArray(ff0Global);
	delete []ffeGlobal;
	
	fffGlobal = new double[N];	
	fffGlobal = generateGlobal(fff);
	MPI_Bcast(fffGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	fff = duplicateArray(fffGlobal);
	delete []fffGlobal;
	
	ffgGlobal = new double[N];	
	ffgGlobal = generateGlobal(ffg);
	MPI_Bcast(ffgGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	ffg = duplicateArray(ffgGlobal);
	delete []ffgGlobal;
	
	ffhGlobal = new double[N];	
	ffhGlobal = generateGlobal(ffh);
	MPI_Bcast(ffhGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	ffh = duplicateArray(ffhGlobal);
	delete []ffhGlobal;
	
	ffiGlobal = new double[N];	
	ffiGlobal = generateGlobal(ffi);
	MPI_Bcast(ffiGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	ffi = duplicateArray(ffiGlobal);
	delete []ffiGlobal;
	
	ffjGlobal = new double[N];	
	ffjGlobal = generateGlobal(ffj);
	MPI_Bcast(ffjGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	ffj = duplicateArray(ffjGlobal);
	delete []ffjGlobal;
	
	ffkGlobal = new double[N];	
	ffkGlobal = generateGlobal(ffk);
	MPI_Bcast(ffkGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	ffk = duplicateArray(ffkGlobal);
	delete []ffkGlobal;
	
	fflGlobal = new double[N];	
	fflGlobal = generateGlobal(ffl);
	MPI_Bcast(fflGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	ffl = duplicateArray(fflGlobal);
	delete []fflGlobal;
	
	
	
	
	gg0Global = new double[N];	
	gg0Global = generateGlobal(gg0);
	MPI_Bcast(gg0Global,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	gg0 = duplicateArray(gg0Global);
	delete []gg0Global;
	
	gg1Global = new double[N];	
	gg1Global = generateGlobal(gg1);
	MPI_Bcast(gg1Global,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	gg1 = duplicateArray(gg1Global);
	delete []gg1Global;
	
	gg2Global = new double[N];	
	gg2Global = generateGlobal(gg2);
	MPI_Bcast(gg2Global,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	gg2 = duplicateArray(gg2Global);
	delete []gg2Global;
	
	gg3Global = new double[N];	
	gg3Global = generateGlobal(gg3);
	MPI_Bcast(gg3Global,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	gg3 = duplicateArray(gg3Global);
	delete []gg3Global;
	
	gg4Global = new double[N];	
	gg4Global = generateGlobal(gg4);
	MPI_Bcast(gg4Global,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	gg4 = duplicateArray(gg4Global);
	delete []gg4Global;
	
	gg5Global = new double[N];	
	gg5Global = generateGlobal(gg5);
	MPI_Bcast(gg5Global,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	gg5 = duplicateArray(gg5Global);
	delete []gg5Global;
	
	gg6Global = new double[N];	
	gg6Global = generateGlobal(gg6);
	MPI_Bcast(gg6Global,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	gg6 = duplicateArray(gg6Global);
	delete []gg6Global;
	
	ggaGlobal = new double[N];	
	ggaGlobal = generateGlobal(gga);
	MPI_Bcast(ggaGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	gga = duplicateArray(ggaGlobal);
	delete []ggaGlobal;
	
	ggbGlobal = new double[N];	
	ggbGlobal = generateGlobal(ggb);
	MPI_Bcast(ggbGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	ggb = duplicateArray(ggbGlobal);
	delete []ggbGlobal;
	
	ggcGlobal = new double[N];	
	ggcGlobal = generateGlobal(ggc);
	MPI_Bcast(ggcGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	ggc = duplicateArray(ggcGlobal);
	delete []ggcGlobal;
	
	ggdGlobal = new double[N];	
	ggdGlobal = generateGlobal(ggd);
	MPI_Bcast(ggdGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	ggd = duplicateArray(ggdGlobal);
	delete []ggdGlobal;
	
	ggeGlobal = new double[N];	
	ggeGlobal = generateGlobal(gge);
	MPI_Bcast(ggeGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	gge = duplicateArray(gg0Global);
	delete []ggeGlobal;
	
	ggfGlobal = new double[N];	
	ggfGlobal = generateGlobal(ggf);
	MPI_Bcast(ggfGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	ggf = duplicateArray(ggfGlobal);
	delete []ggfGlobal;
	
	gggGlobal = new double[N];	
	gggGlobal = generateGlobal(ggg);
	MPI_Bcast(gggGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	ggg = duplicateArray(gggGlobal);
	delete []gggGlobal;
	
	gghGlobal = new double[N];	
	gghGlobal = generateGlobal(ggh);
	MPI_Bcast(gghGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	ggh = duplicateArray(gghGlobal);
	delete []gghGlobal;
	
	ggiGlobal = new double[N];	
	ggiGlobal = generateGlobal(ggi);
	MPI_Bcast(ggiGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	ggi = duplicateArray(ggiGlobal);
	delete []ggiGlobal;
	
	ggjGlobal = new double[N];	
	ggjGlobal = generateGlobal(ggj);
	MPI_Bcast(ggjGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	ggj = duplicateArray(ggjGlobal);
	delete []ggjGlobal;
	
	ggkGlobal = new double[N];	
	ggkGlobal = generateGlobal(ggk);
	MPI_Bcast(ggkGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	ggk = duplicateArray(ggkGlobal);
	delete []ggkGlobal;
	
	gglGlobal = new double[N];	
	gglGlobal = generateGlobal(ggl);
	MPI_Bcast(gglGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	ggl = duplicateArray(gglGlobal);
	delete []gglGlobal;
	
	*/
	
	
	ff0Global = new double[N]; 
	ff1Global = new double[N]; ff2Global = new double[N]; ff3Global = new double[N]; 
	ff4Global = new double[N]; ff5Global = new double[N]; ff6Global = new double[N]; 
	ffaGlobal = new double[N]; ffbGlobal = new double[N]; ffcGlobal = new double[N]; 
	ffdGlobal = new double[N]; ffeGlobal = new double[N]; fffGlobal = new double[N];
	ffgGlobal = new double[N]; ffhGlobal = new double[N]; ffiGlobal = new double[N];
	ffjGlobal = new double[N]; ffkGlobal = new double[N]; fflGlobal = new double[N];
	
	gg0Global = new double[N]; 
	gg1Global = new double[N]; gg2Global = new double[N]; gg3Global = new double[N]; 
	gg4Global = new double[N]; gg5Global = new double[N]; gg6Global = new double[N]; 
	ggaGlobal = new double[N]; ggbGlobal = new double[N]; ggcGlobal = new double[N]; 
	ggdGlobal = new double[N]; ggeGlobal = new double[N]; ggfGlobal = new double[N];
	gggGlobal = new double[N]; gghGlobal = new double[N]; ggiGlobal = new double[N];
	ggjGlobal = new double[N]; ggkGlobal = new double[N]; gglGlobal = new double[N];
	
	// Copy original arrays from different processors to global arrays on ROOT ( arrayGlobal = generateGlobal(array) ) 
	// and distribute them to all processors ( MPI_Bcast(...) )
	
	generateGlobalMask();  
	
	MPI_Bcast(maskGlobal,N,MPI_INT,ROOT,MPI_COMM_WORLD);
	
	
	ff0Global = generateGlobal(ff0);
	ff1Global = generateGlobal(ff1); 
	ff2Global = generateGlobal(ff2);
	ff3Global = generateGlobal(ff3); 
	ff4Global = generateGlobal(ff4);
	ff5Global = generateGlobal(ff5);
	ff6Global = generateGlobal(ff6); 
	ffaGlobal = generateGlobal(ffa);
	ffbGlobal = generateGlobal(ffb); 
	ffcGlobal = generateGlobal(ffc); 
	ffdGlobal = generateGlobal(ffd);
	ffeGlobal = generateGlobal(ffe); 
	fffGlobal = generateGlobal(fff);
	ffgGlobal = generateGlobal(ffg); 
	ffhGlobal = generateGlobal(ffh); 
	ffiGlobal = generateGlobal(ffi);
	ffjGlobal = generateGlobal(ffj); 
	ffkGlobal = generateGlobal(ffk); 
	fflGlobal = generateGlobal(ffl);
	
	gg0Global = generateGlobal(gg0);
	gg1Global = generateGlobal(gg1); 
	gg2Global = generateGlobal(gg2);
	gg3Global = generateGlobal(gg3); 
	gg4Global = generateGlobal(gg4);
	gg5Global = generateGlobal(gg5);
	gg6Global = generateGlobal(gg6); 
	ggaGlobal = generateGlobal(gga);
	ggbGlobal = generateGlobal(ggb); 
	ggcGlobal = generateGlobal(ggc); 
	ggdGlobal = generateGlobal(ggd);
	ggeGlobal = generateGlobal(gge); 
	ggfGlobal = generateGlobal(ggf);
	gggGlobal = generateGlobal(ggg); 
	gghGlobal = generateGlobal(ggh); 
	ggiGlobal = generateGlobal(ggi);
	ggjGlobal = generateGlobal(ggj); 
	ggkGlobal = generateGlobal(ggk); 
	gglGlobal = generateGlobal(ggl);
	
	
	MPI_Bcast(ff0Global,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	MPI_Bcast(ff1Global,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	MPI_Bcast(ff2Global,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	MPI_Bcast(ff3Global,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	MPI_Bcast(ff4Global,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	MPI_Bcast(ff5Global,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	MPI_Bcast(ff6Global,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	MPI_Bcast(ffaGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	MPI_Bcast(ffbGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	MPI_Bcast(ffcGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	MPI_Bcast(ffdGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	MPI_Bcast(ffeGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	MPI_Bcast(fffGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	MPI_Bcast(ffgGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	MPI_Bcast(ffhGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	MPI_Bcast(ffiGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	MPI_Bcast(ffjGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	MPI_Bcast(ffkGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	MPI_Bcast(fflGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	
	MPI_Bcast(gg0Global,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	MPI_Bcast(gg1Global,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	MPI_Bcast(gg2Global,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	MPI_Bcast(gg3Global,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	MPI_Bcast(gg4Global,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	MPI_Bcast(gg5Global,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	MPI_Bcast(gg6Global,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	MPI_Bcast(ggaGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	MPI_Bcast(ggbGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	MPI_Bcast(ggcGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	MPI_Bcast(ggdGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	MPI_Bcast(ggeGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	MPI_Bcast(ggfGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	MPI_Bcast(gggGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	MPI_Bcast(gghGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	MPI_Bcast(ggiGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	MPI_Bcast(ggjGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	MPI_Bcast(ggkGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	MPI_Bcast(gglGlobal,N,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	
	// Now all old arrays can re-initialised; in initialise LX will be the full value, 
	// so all arrays and LX, N are twice the size from now on
	
	initialise(); 
		
	cout << "Proces " << rank <<": re-initialised" << endl;
	cout << "Proces " << rank << ": new LX = " << LX << endl;
	
	// Copy all ff,gg,mask values from global arrays to new, double-sized arrays;
	// in new half mirror-image of old half
	
	mask = duplicateArray_int(maskGlobal);
	
	ff0 = duplicateArray(ff0Global);
	ff1 = duplicateArray(ff1Global);
	ff2 = duplicateArray(ff2Global);
	ff3 = duplicateArray(ff3Global);
	ff4 = duplicateArray(ff4Global);
	ff5 = duplicateArray(ff5Global);
	ff6 = duplicateArray(ff6Global);
	ffa = duplicateArray(ffaGlobal);
	ffb = duplicateArray(ffbGlobal);
	ffc = duplicateArray(ffcGlobal);
	ffd = duplicateArray(ffdGlobal);
	ffe = duplicateArray(ffeGlobal);
	fff = duplicateArray(fffGlobal);
	ffg = duplicateArray(ffgGlobal);
	ffh = duplicateArray(ffhGlobal);
	ffi = duplicateArray(ffiGlobal);
	ffj = duplicateArray(ffjGlobal);
	ffk = duplicateArray(ffkGlobal);
	ffl = duplicateArray(fflGlobal);
	
	gg0 = duplicateArray(gg0Global);
	gg1 = duplicateArray(gg1Global);
	gg2 = duplicateArray(gg2Global);
	gg3 = duplicateArray(gg3Global);
	gg4 = duplicateArray(gg4Global);
	gg5 = duplicateArray(gg5Global);
	gg6 = duplicateArray(gg6Global);
	gga = duplicateArray(ggaGlobal);
	ggb = duplicateArray(ggbGlobal);
	ggc = duplicateArray(ggcGlobal);
	ggd = duplicateArray(ggdGlobal);
	gge = duplicateArray(ggeGlobal);
	ggf = duplicateArray(ggfGlobal);
	ggg = duplicateArray(gggGlobal);
	ggh = duplicateArray(gghGlobal);
	ggi = duplicateArray(ggiGlobal);
	ggj = duplicateArray(ggjGlobal);
	ggk = duplicateArray(ggkGlobal);
	ggl = duplicateArray(gglGlobal);
	
	// All done, can de-allocate temporary arrays
		
	delete []ff0Global;
	delete []ff1Global;
	delete []ff2Global;
	delete []ff3Global;
	delete []ff4Global;
	delete []ff5Global;
	delete []ff6Global;
	delete []ffaGlobal;
	delete []ffbGlobal;
	delete []ffcGlobal;
	delete []ffdGlobal;
	delete []ffeGlobal;
	delete []fffGlobal;
	delete []ffgGlobal;
	delete []ffhGlobal;
	delete []ffiGlobal;
	delete []ffjGlobal;
	delete []ffkGlobal;
	delete []fflGlobal;
	
	delete []gg0Global;
	delete []gg1Global;
	delete []gg2Global;
	delete []gg3Global;
	delete []gg4Global;
	delete []gg5Global;
	delete []gg6Global;
	delete []ggaGlobal;
	delete []ggbGlobal;
	delete []ggcGlobal;
	delete []ggdGlobal;
	delete []ggeGlobal;
	delete []ggfGlobal;
	delete []gggGlobal;
	delete []gghGlobal;
	delete []ggiGlobal;
	delete []ggjGlobal;
	delete []ggkGlobal;
	delete []gglGlobal;
	
	
	delete []maskGlobal;
	delete []pGlobal;	//still smaller size
	//delete []nGlobal;
	
	maskGlobal = new int[N];
	pGlobal = new double[N]; //correct post-equil size
	//nGlobal = new double[N];
	
	for (k=k1; k< k2; k++)
		relabel();
	for (k=k1; k< k2; k++)
		leakSearch();
	for (k=k1; k< k2; k++){
		if (mask[k] == 28) {n[k] = 1.0; p[k] = 0.0;}
	}
	
	generateGlobalMask();
	exchangeMask();
	
	/*
	if(rank==ROOT)
	{
		int i, j, h;	
		
		char filename[25];
		sprintf(filename, "./densitydata/gg0_dupl%ld.m", t);			//Create a name for file that contain data 
		//sprintf(filename, "./dd%ld.m", t);
		
		ofstream file(filename);
		file.precision(4);
		for( h = 0 ; h < LZ ; h++) 
		{   
			file << "gg0" << t << "(:,:," << h+1 << ")=[" << endl;
			for( i = 0 ; i < LX/2 ; i++) 
			{
				for( j = 0 ; j < LY ; j++) 
				{
					k = h + j*LZ + i*LY*LZ;
					//if( maskGlobal[k]==28) 
					//	file << -2.0 << " ";  
					//else
						file << gg0Global[k] << " ";
					
				}
				file << endl;
			}
			file <<"];" << endl;
		}
		file.close();
	}
	
	*/		
	
}
	
