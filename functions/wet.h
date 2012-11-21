//Marco, 29 june 2011, Lisa, 28 July 2011
//Costruction of a class for Lattice Boltzmann wetting simulation

#include <cstring>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <mpi.h>
#include "String.hh"
#include <string>
#include <sstream>
//#include "RT_Timer.hh"

using namespace std;

#define NRATIOFLUID 1000	//number of relaxation times stored considering ratio of density of two fluid
#define ROOT 0				//definig root process
const int M[19][19]={{  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1},
			     {-30,-11,-11,-11,-11,-11,-11,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8},
			     { 12, -4, -4, -4, -4, -4, -4,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1},
			     {  0,  1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  0,  0,  0,  0,  1, -1,  1, -1},
			     {  0, -4,  4,  0,  0,  0,  0,  1, -1,  1, -1,  0,  0,  0,  0,  1, -1,  1, -1},
			     {  0,  0,  0,  1, -1,  0,  0,  1,  1, -1, -1,  1, -1,  1, -1,  0,  0,  0,  0},
			     {  0,  0,  0, -4,  4,  0,  0,  1,  1, -1, -1,  1, -1,  1, -1,  0,  0,  0,  0},
			     {  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1,  1, -1, -1,  1,  1, -1, -1},
			     {  0,  0,  0,  0,  0, -4,  4,  0,  0,  0,  0,  1,  1, -1, -1,  1,  1, -1, -1},
			     {  0,  2,  2, -1, -1, -1, -1,  1,  1,  1,  1, -2, -2, -2, -2,  1,  1,  1,  1},
			     {  0, -4, -4,  2,  2,  2,  2,  1,  1,  1,  1, -2, -2, -2, -2,  1,  1,  1,  1},
			     {  0,  0,  0,  1,  1, -1, -1,  1,  1,  1,  1,  0,  0,  0,  0, -1, -1, -1, -1},
			     {  0,  0,  0, -2, -2,  2,  2,  1,  1,  1,  1,  0,  0,  0,  0, -1, -1, -1, -1},
			     {  0,  0,  0,  0,  0,  0,  0,  1, -1, -1,  1,  0,  0,  0,  0,  0,  0,  0,  0},
			     {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1, -1, -1,  1,  0,  0,  0,  0},
			     {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1, -1, -1,  1},
			     {  0,  0,  0,  0,  0,  0,  0,  1, -1,  1, -1,  0,  0,  0,  0, -1,  1, -1,  1},
			     {  0,  0,  0,  0,  0,  0,  0, -1, -1,  1,  1,  1, -1,  1, -1,  0,  0,  0,  0},
			     {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1, -1,  1,  1,  1,  1, -1, -1}};

class wet
{
	int size;

	long t; 				// index 
	
	//bool equilfirst;		//two drops in the first place, or one drop that is equilibrated and then duplicated, settings.par
	bool velocityprofile;
	
	
	int dimensions;			//quasi-2d or 3d, defined in wet.par
	int LX, LY, LZ, LXc;			//Number of lattice points, defined in wet.par
	
	//int nbEqStep;			//Number of total iterations, define in wet.par	
 	//int equilTime;			//time of equilibrium, defined in wet.par				(SISTEMA)
 	int writeStep;  		//Number of steps at which output is written, defined in wet.par
	int infoStep; 			//Number of steps at witch info is collected, defined in wet.par
	double teta1, teta2, tetac;    //, teta_CB_fix;		//Contact angles, defined in wet.par	
	double tauliquid, taugas;	//Relaxation times of liquid and gas, defined in wet.par
	double tau0, tau1;		//Relaxation time 	
	double tausurf;		    //Relaxation time for fixing phi-derivative at surface, def. in wet.par, used in densitiesAtSolidBoundaries() 
	double kappa_p;
	double A, B, gama; 		//A, Gama (diffusion) defined in wet.par
	double G[3];			//forces acting on lattice points, defined in wet.par used in initialise()  
	double p_thresh;        //threshold in order parameter for differenciation interface <-> surface <-> bulk (not really used, I mostly look at 'natural threshold' = tanh(1))

	
	char geometry[15];  //Substrate geometry, defined in substrate.par
	int Dx, Dy;			//Dimension of a reactangular post (in lattice points), defined in substrate.par
 	int Dh;				//Height of posts or height of floor (in lattice points),  defined in substrate.par
    int Px , Py ;       // X and Y Position of thread
 	int PeriodX, PeriodY;			//Period of posts (<Dx and <Dy) (in lattice points), defined in substrate.par
	int PostSymm; 		//Posts with gap at LX/2 (0) or with post (1)?, defined in substrate.par and used in initialiseSurface()
	

	double dropletR;			//Radius of a free droplet, defined in spherical.par
	double dropletCenterX;		//X position of free droplet (in lattice points), defined in spherical.par
	double dropletCenterY;		//Y position of free droplet (in lattice points), defined in spherical.par	
	double dropletCenterZ;		//Z position of free droplet (in lattice points), defined in spherical.par
	double initUX, initUY, initUZ;	//Initial velocity of the droplet, defined in spherical.par used in LGConfig()

	double teta_CB, areaFrac; 
	
	double dx;			//Lattice step (initialised in initialise() dx = 1.0), used in equilibrium(). To change this must change derivatives
	double dt;			//Lattice time step (initialised in initialise() dt = 1.0), in equilibrium(). To change this must change derivatives
	double c, c2;		//Evolution velocity (dx/dt) and its square (initialised in initialise()), used in equilibrium()

	double z1, z2, z3, z4, z5, z6, z7;	//constant used to calculate equilibrium function, initialise in initialise() and used in equilibrium() 	

	int N;					//Number of total node (LX*LY*LZ), in initialise().
	long k, k1, k2;			//Index of current lattice point (between 1 and N=LX*LY*LZ)
	int xk, yk, zk;			//lattice points, 3-D index version of index k, in ComputeCoordinate()

	int processN;			//Number of total nodes per process, in initialise();

	int leftProcess, rightProcess;	//identification numbers for neighboor processes in the ring
			
	
	//ff(i) is a pointer to a double vector of dimension N for particle of fluid f, with direction i, after streaming and collision
	//They are initialize in initialize()
	double *ff0, *ff1, *ff2, *ff3, *ff4, *ff5, *ff6, *ffa, *ffb, *ffc, *ffd, *ffe, *fff, *ffg, *ffh, *ffi, *ffj, *ffk, *ffl; 
	//gg(i) is a pointer to a double vector of dimension N for particle of fluid g, with direction i, after streaming and collision
	//They are initialize in initialize()
	double *gg0, *gg1, *gg2, *gg3, *gg4, *gg5, *gg6, *gga, *ggb, *ggc, *ggd, *gge, *ggf, *ggg, *ggh, *ggi, *ggj, *ggk, *ggl; 
	//fn(i) is a pointer to a double vector of dimension N for particle of fluid f, with direction i, after collision, before streaming
	//They are initialize in initialize()
	double *fn0, *fn1, *fn2, *fn3, *fn4, *fn5, *fn6, *fna, *fnb, *fnc, *fnd, *fne, *fnf, *fng, *fnh, *fni, *fnj, *fnk, *fnl; 
	//gn(i) is a pointer to a double vector of dimension N for particle of fluid f, with direction i, after collision, before streaming
	//They are initialize in initialize();
	double *gn0, *gn1, *gn2, *gn3, *gn4, *gn5, *gn6, *gna, *gnb, *gnc, *gnd, *gne, *gnf, *gng, *gnh, *gni, *gnj, *gnk, *gnl; 

	//Pointers to a N dimensional vector that represents TOTAL density of fluid f (*n) an g (*p) at lattice point k,
	//created in initialise(), initialised in LGConfig() and calculated in computeMomenta()
	double *n, *p, *pGlobal, *nGlobal; 
    //double *ndup, *pdup;
	//Pointers to a vector of dimendion N that represent velocity of fluid at lattice point, in initialise() and LGConfig()	      (SISTEMA)	
	
	//double *uxs, *uys, *uzs;
	double *uxs, *uys, *uzs, *uGlobal;		
	//Pointers to a vector of dimendion N that represent forces acting on lattice points, in initialise() and collision()	      (SISTEMA)
	double  *forcex, *forcey, *forcez;
	//Pointers to a vector of dimension N that is the neighbours of each lattice point, in initialise()
	long *dd1, *dd2, *dd3, *dd4, *dd5, *dd6, *dda, *ddb, *ddc, *ddd, *dde, *ddf, *ddg, *ddh, *ddi, *ddj, *ddk, *ddl;       
	//Pointers to a vector of dimension N that is lattice ID for topolocially patterned substrate, defined in initialise()
	//usually mask=0 for bulk fluid, mask=1 for fluid near solid surface, mask=28 for surface (that are arbitrary parameters)
	int *mask, *maskGlobal;  
	
	double *dissipation, *dissGlobal;   //array of dissipation values for dissipation from shear, calculated in computeDissipation()
	double *diffDiss, *diffDissGlobal;  //array of dissipation values for dissipation from diffusion, calculated in computeDissipation()
	double *mu, *muGlobal ;                         //array for chemical potential, needed for dissipation from diffusion, calc. in equilibrium(), used in computeDissipation()
	//double *eBal, *eBal_old, *eBalDiff, *eBalDiffGlobal;    //FOR NOW JUST TESTING: dissipation from energy balance at every lattice node   
	
	double dp_dx, dp_dy, dp_dz;	//Density gradient of order parameter, used in initialise(), calculated in collision()
	double del2n, del2p;		//Laplacian of density, used in equilibrium(), calculated in collision()
	double ux, uy, uz;		//Temporary velocity variables, in equilibrium() and collision()						
	double ux2, uy2, uz2;		//Square of temporary velocity variable, and total square velocity, in equilibrium()
	double uxuy2, uyuz2, uxuz2; 	//variables used to calculate equilibrium functions, in equilibrium()
		
	double fbulk1, fbulk2;		//variables used to calculate equilibrium function of f, in equilibrium()
	double nz1, nz2, nz5;		//variables used to calculate equilibrium function of f, in equilibrium()
	double nz6ux2, nz6uy2, nz6uz2;	//variables used to calculate equilibrium function of f, in equilibrium()
	double nz7uxuy, nz7uyuz, nz7uxuz;	//variables used to calculate equilibrium function of f, in equilibrium()
	double a1, a2, a3, a4, a5, a6;	//variables used to calculate equilibrium function of f, in equilibrium()
	double a12, a23, a13;		//variables used to calculate equilibrium function of f, in equilibrium()

	double gbulk1, gbulk2;      	//variables used to calculate equilibrium function of g, in equilibrium()                  
	double pz1, pz2, pz5;		//variables used to calculate equilibrium function of g, in equilibrium()	
	double pz6ux2, pz6uy2, pz6uz2;	//variables used to calculate equilibrium function of g, in equilibrium()
	double pz7uxuy, pz7uyuz, pz7uxuz;	//variables used to calculate equilibrium function of g, in equilibrium()
	double fx, fy, fz;		//variables used to calculate force on fluid, in equilibrium()
	double Mxxyy, Myyzz, Mxxzz;	//variables used to calculate force on fluid, in equilibrium()
	double Mxyyx, Myzzy, Mxzzx;	//variables used to calculate force on fluid, in equilibrium()
	double fsq;			//variable used to calculate force on fluid, in equilibrium()
	double tau2;			//relaxation time between fluid 1 and 2 used to calculate force on fluid, in equilibrium()	
	    
	//fe(i) is the equilibrium partial density in direction i of fluid f, in equilibrium()   
	double fe0, fe1, fe2, fe3, fe4, fe5, fe6, fea, feb, fec, fed, fee, fef, feg, feh, fei, fej, fek, fel; 
	//ge(i) is the equilibrium partial density in direction i of fluid g, in equilibrium()
	double ge0, ge1, ge2, ge3, ge4, ge5, ge6, gea, geb, gec, ged, gee, gef, geg, geh, gei, gej, gek, gel;

	
	double phi11,phi12;		//Surface free energy parameters, calculated in initialise() (see literature)

	double energy, energy_n, energy_init, energy_n_init;					//energy, used in LBAlgorithm()						
	double bulkE, bulkE_n, bulkE_init, bulkE_n_init;						//bulk energy, used in LBAlgorithm()					
	double interfaceE, interfaceE_n, interfaceE_init, interfaceE_n_init;	//interface energy, used in LBAlgorithm()				
	double surfaceE, surfaceE_n, surfaceE_init, surfaceE_n_init;		    //surface energy, used in LBAlgorithm()	
	double p_thresh_n;														//threshold to distinguish between bulk and 
	double surfArea, surfArea_init, LFree, LFree_init;
	int detachflag;
	double clLength, clLength_n, clLength_init, clLength_n_init;
	
	double dissSum_liq, dissSum_gas, dissSum_int;
	double dissInt_liq, dissInt_gas, dissInt_int;
	double dissInt_liq_old, dissInt_gas_old, dissInt_int_old;  //dissipation()
	double diffDissSum_liq, diffDissSum_gas, diffDissSum_int, diffDiss_init;
	double diffDissInt_liq, diffDissInt_gas, diffDissInt_int;
	double diffDissInt_liq_old, diffDissInt_gas_old, diffDissInt_int_old;  //dissipation()
	double mLiq, mGas, mInt;  //amount of liquid, interface, gas
	
	double energy_n_old, kin_old; //dissipation()
	
    double  Oh,we,w,Re,vg,mo,stg,vl,stl; //Quantities used in Dimensionless num
	
	double RXcm, RXcmOld;		//Center of mass of fluid f and the same "old" variable along x, in ARolling()
	double RZcm;				//Center of mass of fluid f along z, in ARolling()
	double RNtot;				//TOTAL mass, sum of n for all N lattice points, calculated in ARolling()
	double RVXcm, RVYcm, RVZcm;	//Velocity of center of mass of fluid f, in ARolling()
	double kin, kin_init, kinG, kinL;					//TOTAL kinetic energy of fluid f, calculated in ARolling()
	double kin_CM;				//kinetic energy of center of mass of fluid f, calculated in ARolling()	
	double ratio_kin;			//ratio between kinetic energy of center of mass and total kinetic energy (kin_CM/kin), in ARolling()
	double eccentricity;
	double d_xk, d_yk, d_zk;	//lattice points, 3-D index version of index k (like xk, yk, zk, but double), in ARolling()
	double RNtot_gas;
	
	//bool afterequilflag;
	
	double  force0, force1, force2, force3, force4, force5, force6, forcea, forceb, forcec, forced, forcee, forcef, forceg, forceh, forcei, forcej, forcek, forcel;
	//double gradsq, en;		//square of density gradient (grad(p)*grad(p)), and used in computeFreeEnergy()
	
	double M00[NRATIOFLUID], M01[NRATIOFLUID], M02[NRATIOFLUID], M03[NRATIOFLUID], M04[NRATIOFLUID], M05[NRATIOFLUID], M06[NRATIOFLUID];
	double M0a[NRATIOFLUID], M0b[NRATIOFLUID], M0c[NRATIOFLUID], M0d[NRATIOFLUID], M0e[NRATIOFLUID], M0f[NRATIOFLUID], M0g[NRATIOFLUID]; 
	double M0h[NRATIOFLUID], M0i[NRATIOFLUID], M0j[NRATIOFLUID], M0k[NRATIOFLUID], M0l[NRATIOFLUID];
	double M10[NRATIOFLUID], M11[NRATIOFLUID], M12[NRATIOFLUID], M13[NRATIOFLUID], M14[NRATIOFLUID], M15[NRATIOFLUID], M16[NRATIOFLUID];
	double M1a[NRATIOFLUID], M1b[NRATIOFLUID], M1c[NRATIOFLUID], M1d[NRATIOFLUID], M1e[NRATIOFLUID], M1f[NRATIOFLUID], M1g[NRATIOFLUID]; 
	double M1h[NRATIOFLUID], M1i[NRATIOFLUID], M1j[NRATIOFLUID], M1k[NRATIOFLUID], M1l[NRATIOFLUID];
	double M20[NRATIOFLUID], M21[NRATIOFLUID], M22[NRATIOFLUID], M23[NRATIOFLUID], M24[NRATIOFLUID], M25[NRATIOFLUID], M26[NRATIOFLUID];
	double M2a[NRATIOFLUID], M2b[NRATIOFLUID], M2c[NRATIOFLUID], M2d[NRATIOFLUID], M2e[NRATIOFLUID], M2f[NRATIOFLUID], M2g[NRATIOFLUID]; 
	double M2h[NRATIOFLUID], M2i[NRATIOFLUID], M2j[NRATIOFLUID], M2k[NRATIOFLUID], M2l[NRATIOFLUID];
	double M30[NRATIOFLUID], M31[NRATIOFLUID], M32[NRATIOFLUID], M33[NRATIOFLUID], M34[NRATIOFLUID], M35[NRATIOFLUID], M36[NRATIOFLUID];
	double M3a[NRATIOFLUID], M3b[NRATIOFLUID], M3c[NRATIOFLUID], M3d[NRATIOFLUID], M3e[NRATIOFLUID], M3f[NRATIOFLUID], M3g[NRATIOFLUID];
	double M3h[NRATIOFLUID], M3i[NRATIOFLUID], M3j[NRATIOFLUID], M3k[NRATIOFLUID], M3l[NRATIOFLUID];
	double M40[NRATIOFLUID], M41[NRATIOFLUID], M42[NRATIOFLUID], M43[NRATIOFLUID], M44[NRATIOFLUID], M45[NRATIOFLUID], M46[NRATIOFLUID];
	double M4a[NRATIOFLUID], M4b[NRATIOFLUID], M4c[NRATIOFLUID], M4d[NRATIOFLUID], M4e[NRATIOFLUID], M4f[NRATIOFLUID], M4g[NRATIOFLUID];
	double M4h[NRATIOFLUID], M4i[NRATIOFLUID], M4j[NRATIOFLUID], M4k[NRATIOFLUID], M4l[NRATIOFLUID];
	double M50[NRATIOFLUID], M51[NRATIOFLUID], M52[NRATIOFLUID], M53[NRATIOFLUID], M54[NRATIOFLUID], M55[NRATIOFLUID], M56[NRATIOFLUID];
	double M5a[NRATIOFLUID], M5b[NRATIOFLUID], M5c[NRATIOFLUID], M5d[NRATIOFLUID], M5e[NRATIOFLUID], M5f[NRATIOFLUID], M5g[NRATIOFLUID]; 
	double M5h[NRATIOFLUID], M5i[NRATIOFLUID], M5j[NRATIOFLUID], M5k[NRATIOFLUID], M5l[NRATIOFLUID];
	double M60[NRATIOFLUID], M61[NRATIOFLUID], M62[NRATIOFLUID], M63[NRATIOFLUID], M64[NRATIOFLUID], M65[NRATIOFLUID], M66[NRATIOFLUID];
	double M6a[NRATIOFLUID], M6b[NRATIOFLUID], M6c[NRATIOFLUID], M6d[NRATIOFLUID], M6e[NRATIOFLUID], M6f[NRATIOFLUID], M6g[NRATIOFLUID];
	double M6h[NRATIOFLUID], M6i[NRATIOFLUID], M6j[NRATIOFLUID], M6k[NRATIOFLUID], M6l[NRATIOFLUID];
	double Ma0[NRATIOFLUID], Ma1[NRATIOFLUID], Ma2[NRATIOFLUID], Ma3[NRATIOFLUID], Ma4[NRATIOFLUID], Ma5[NRATIOFLUID], Ma6[NRATIOFLUID];
	double Maa[NRATIOFLUID], Mab[NRATIOFLUID], Mac[NRATIOFLUID], Mad[NRATIOFLUID], Mae[NRATIOFLUID], Maf[NRATIOFLUID], Mag[NRATIOFLUID];
	double Mah[NRATIOFLUID], Mai[NRATIOFLUID], Maj[NRATIOFLUID], Mak[NRATIOFLUID], Mal[NRATIOFLUID];
	double Mb0[NRATIOFLUID], Mb1[NRATIOFLUID], Mb2[NRATIOFLUID], Mb3[NRATIOFLUID], Mb4[NRATIOFLUID], Mb5[NRATIOFLUID], Mb6[NRATIOFLUID];
	double Mba[NRATIOFLUID], Mbb[NRATIOFLUID], Mbc[NRATIOFLUID], Mbd[NRATIOFLUID], Mbe[NRATIOFLUID], Mbf[NRATIOFLUID], Mbg[NRATIOFLUID];
	double Mbh[NRATIOFLUID], Mbi[NRATIOFLUID], Mbj[NRATIOFLUID], Mbk[NRATIOFLUID], Mbl[NRATIOFLUID];
	double Mc0[NRATIOFLUID], Mc1[NRATIOFLUID], Mc2[NRATIOFLUID], Mc3[NRATIOFLUID], Mc4[NRATIOFLUID], Mc5[NRATIOFLUID], Mc6[NRATIOFLUID];
	double Mca[NRATIOFLUID], Mcb[NRATIOFLUID], Mcc[NRATIOFLUID], Mcd[NRATIOFLUID], Mce[NRATIOFLUID], Mcf[NRATIOFLUID], Mcg[NRATIOFLUID];
	double Mch[NRATIOFLUID], Mci[NRATIOFLUID], Mcj[NRATIOFLUID], Mck[NRATIOFLUID], Mcl[NRATIOFLUID];
	double Md0[NRATIOFLUID], Md1[NRATIOFLUID], Md2[NRATIOFLUID], Md3[NRATIOFLUID], Md4[NRATIOFLUID], Md5[NRATIOFLUID], Md6[NRATIOFLUID];
	double Mda[NRATIOFLUID], Mdb[NRATIOFLUID], Mdc[NRATIOFLUID], Mdd[NRATIOFLUID], Mde[NRATIOFLUID], Mdf[NRATIOFLUID], Mdg[NRATIOFLUID];
	double Mdh[NRATIOFLUID], Mdi[NRATIOFLUID], Mdj[NRATIOFLUID], Mdk[NRATIOFLUID], Mdl[NRATIOFLUID];
	double Me0[NRATIOFLUID], Me1[NRATIOFLUID], Me2[NRATIOFLUID], Me3[NRATIOFLUID], Me4[NRATIOFLUID], Me5[NRATIOFLUID], Me6[NRATIOFLUID];
	double Mea[NRATIOFLUID], Meb[NRATIOFLUID], Mec[NRATIOFLUID], Med[NRATIOFLUID], Mee[NRATIOFLUID], Mef[NRATIOFLUID], Meg[NRATIOFLUID];
	double Meh[NRATIOFLUID], Mei[NRATIOFLUID], Mej[NRATIOFLUID], Mek[NRATIOFLUID], Mel[NRATIOFLUID];
	double Mf0[NRATIOFLUID], Mf1[NRATIOFLUID], Mf2[NRATIOFLUID], Mf3[NRATIOFLUID], Mf4[NRATIOFLUID], Mf5[NRATIOFLUID], Mf6[NRATIOFLUID];
	double Mfa[NRATIOFLUID], Mfb[NRATIOFLUID], Mfc[NRATIOFLUID], Mfd[NRATIOFLUID], Mfe[NRATIOFLUID], Mff[NRATIOFLUID], Mfg[NRATIOFLUID];
	double Mfh[NRATIOFLUID], Mfi[NRATIOFLUID], Mfj[NRATIOFLUID], Mfk[NRATIOFLUID], Mfl[NRATIOFLUID];
	double Mg0[NRATIOFLUID], Mg1[NRATIOFLUID], Mg2[NRATIOFLUID], Mg3[NRATIOFLUID], Mg4[NRATIOFLUID], Mg5[NRATIOFLUID], Mg6[NRATIOFLUID];
	double Mga[NRATIOFLUID], Mgb[NRATIOFLUID], Mgc[NRATIOFLUID], Mgd[NRATIOFLUID], Mge[NRATIOFLUID], Mgf[NRATIOFLUID], Mgg[NRATIOFLUID];
	double Mgh[NRATIOFLUID], Mgi[NRATIOFLUID], Mgj[NRATIOFLUID], Mgk[NRATIOFLUID], Mgl[NRATIOFLUID];
	double Mh0[NRATIOFLUID], Mh1[NRATIOFLUID], Mh2[NRATIOFLUID], Mh3[NRATIOFLUID], Mh4[NRATIOFLUID], Mh5[NRATIOFLUID], Mh6[NRATIOFLUID];
	double Mha[NRATIOFLUID], Mhb[NRATIOFLUID], Mhc[NRATIOFLUID], Mhd[NRATIOFLUID], Mhe[NRATIOFLUID], Mhf[NRATIOFLUID], Mhg[NRATIOFLUID];
	double Mhh[NRATIOFLUID], Mhi[NRATIOFLUID], Mhj[NRATIOFLUID], Mhk[NRATIOFLUID], Mhl[NRATIOFLUID];
	double Mi0[NRATIOFLUID], Mi1[NRATIOFLUID], Mi2[NRATIOFLUID], Mi3[NRATIOFLUID], Mi4[NRATIOFLUID], Mi5[NRATIOFLUID], Mi6[NRATIOFLUID];
	double Mia[NRATIOFLUID], Mib[NRATIOFLUID], Mic[NRATIOFLUID], Mid[NRATIOFLUID], Mie[NRATIOFLUID], Mif[NRATIOFLUID], Mig[NRATIOFLUID];
	double Mih[NRATIOFLUID], Mii[NRATIOFLUID], Mij[NRATIOFLUID], Mik[NRATIOFLUID], Mil[NRATIOFLUID];
	double Mj0[NRATIOFLUID], Mj1[NRATIOFLUID], Mj2[NRATIOFLUID], Mj3[NRATIOFLUID], Mj4[NRATIOFLUID], Mj5[NRATIOFLUID], Mj6[NRATIOFLUID];
	double Mja[NRATIOFLUID], Mjb[NRATIOFLUID], Mjc[NRATIOFLUID], Mjd[NRATIOFLUID], Mje[NRATIOFLUID], Mjf[NRATIOFLUID], Mjg[NRATIOFLUID];
	double Mjh[NRATIOFLUID], Mji[NRATIOFLUID], Mjj[NRATIOFLUID], Mjk[NRATIOFLUID], Mjl[NRATIOFLUID];
	double Mk0[NRATIOFLUID], Mk1[NRATIOFLUID], Mk2[NRATIOFLUID], Mk3[NRATIOFLUID], Mk4[NRATIOFLUID], Mk5[NRATIOFLUID], Mk6[NRATIOFLUID];
	double Mka[NRATIOFLUID], Mkb[NRATIOFLUID], Mkc[NRATIOFLUID], Mkd[NRATIOFLUID], Mke[NRATIOFLUID], Mkf[NRATIOFLUID], Mkg[NRATIOFLUID];
	double Mkh[NRATIOFLUID], Mki[NRATIOFLUID], Mkj[NRATIOFLUID], Mkk[NRATIOFLUID], Mkl[NRATIOFLUID];
	double Ml0[NRATIOFLUID], Ml1[NRATIOFLUID], Ml2[NRATIOFLUID], Ml3[NRATIOFLUID], Ml4[NRATIOFLUID], Ml5[NRATIOFLUID], Ml6[NRATIOFLUID];
	double Mla[NRATIOFLUID], Mlb[NRATIOFLUID], Mlc[NRATIOFLUID], Mld[NRATIOFLUID], Mle[NRATIOFLUID], Mlf[NRATIOFLUID], Mlg[NRATIOFLUID];
	double Mlh[NRATIOFLUID], Mli[NRATIOFLUID], Mlj[NRATIOFLUID], Mlk[NRATIOFLUID], Mll[NRATIOFLUID];

	double df0, df1, df2, df3, df4, df5, df6, dfa, dfb, dfc, dfd, dfe, dff, dfg, dfh, dfi, dfj, dfk, dfl;

	int l;				//ratio between densities of two fluids
	

	
	void computeCoordinate(int); 	//computes the 3-D coordinate of index k, used in initialise()
	void initialiseSurface(void);	//Initialise the surface geometries, used in initialise()
	void LGConfig(void);			//Initialise the Binary fluid configuration, used in initialise()
	void relabel(void);				//rename the mask of fluid (mask=0) near a solid surface with mask=1 , used  in initialise()
	void leakSearch(void);			//checks that all fluid near the surface is been counted, used in initialise()
	void equilibrium(void);			//calculates the equilibrium dist function, used in initialise()
	double sign(double value);		//returns the sign of some input value, used in initialise() 
	void computeMomenta(void);		//computes n (density), phi (energy) and the velocity, used in LBAlgorithm()
	void ARolling(void);			//computes dynimic variables (kinetic energy, center of mass velocity,...),  used in LBAlgorithm()
	void densitiesAtSolidBoundaries(void);	//set the density p of fluid g at boundaries
	void computeFreeEnergy(void);
	void invComputeCoordinate(void);	//calculates the k index from the node indeces
	void collision(void);   
	void writeDensityFile(void);
	void writeZPlanXVelocityFile(const int);
	void saveFiles(void);
	void propagation(void);
	void applyBoundaryConditions(void);
	int  ReadInput(void);			//Read data about drop, surface and fluid features
	void initialise(void);			//Initialise computational variables and problem geometry
	void makematrix(void);
	void commReadData(void);
	void exchangeMask(void);
	void exchangePhi(void);
	void generateGlobalMask(void);
	void generatePhiGlobal(void);
	void generateNGlobal(void);
	void exchangeDensities(void);
	double* generateGlobal(double*);
	double* duplicateArray(double*);
	int* duplicateArray_int(int*);
	void computeContactArea(void);
	void writeInfoFile(void);
	void writeJumpFile(void);
	void computeDissipation(void);
	void writeDissipationFile(void);
	void exchangeVelocities(void);
	void exchangeChemPot(void);
	void LGConfigRev(void);
	void exchangeDensities_ffgg(void);
    void Dimensionlessnum(void);
	public:
		
		int rank;
		int equilTime;
		int floorTime;
		int nbEqStep;
		bool afterequilflag;
		bool equilfirst;
		int duplicationtype;
		double teta_CB_fix;
	
	
	
		wet(void);				//Costructor (i.e. initialize variable define in the class)
		void LBAlgorithm(int, int);	//Lattice Boltzmann algorithm
		//void initialise(void);
	    void duplicateDrop1(void);
		void duplicateDrop2(void);
		void initialiseSurfaceLater(void);
		~wet(void);				//Destructor (i.e. free memory when structure is no more used)
};
