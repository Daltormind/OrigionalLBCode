//Destructor: clear memory
//Marco, 28 June 2011

#include "wet.h"
#include <iostream>

using namespace std;

wet::~wet(void)
{	
	
	//Deallocate memory for variables                           
	//ff(i) is a pointer to a double vector of dimension N for particle of fluid f, with direction i, after streaming and collision
	delete []ff0; 
	delete []ff1; delete []ff2; delete []ff3; 
	delete []ff4; delete []ff5; delete []ff6; 
	delete []ffa; delete []ffb; delete []ffc; 
	delete []ffd; delete []ffe; delete []fff;
	delete []ffg; delete []ffh; delete []ffi;
	delete []ffj; delete []ffk; delete []ffl;

	//gg(i) is a pointer to a double vector of dimension N for particle of fluid g, with direction i, after streaming and collision
	delete []gg0; 
	delete []gg1; delete []gg2; delete []gg3; 
	delete []gg4; delete []gg5; delete []gg6; 
	delete []gga; delete []ggb; delete []ggc; 
	delete []ggd; delete []gge; delete []ggf;
	delete []ggg; delete []ggh; delete []ggi;
	delete []ggj; delete []ggk; delete []ggl;

	//fn(i) is a pointer to a double vector of dimension N for particle of fluid f, with direction i, after collision, before streaming
	delete []fn0; 
	delete []fn1; delete []fn2; delete []fn3; 
	delete []fn4; delete []fn5; delete []fn6; 
	delete []fna; delete []fnb; delete []fnc; 
	delete []fnd; delete []fne; delete []fnf;
	delete []fng; delete []fnh; delete []fni;
	delete []fnj; delete []fnk; delete []fnl;
	
	//gn(i) is a pointer to a double vector of dimension N for particle of fluid f, with direction i, after collision, before streaming
	delete []gn0; 
	delete []gn1; delete []gn2; delete []gn3; 
	delete []gn4; delete []gn5; delete []gn6; 
	delete []gna; delete []gnb; delete []gnc; 
	delete []gnd; delete []gne; delete []gnf;
	delete []gng; delete []gnh; delete []gni;
	delete []gnj; delete []gnk; delete []gnl;
  
	delete []n;			//Density of lattice points		(SISTEMA)
	delete []p;			//Density of lattice points		(SISTEMA)

	delete []uxs; 	//Velocity of lattice points			(SISTEMA)
	delete []uys; 
	delete []uzs; 
  	
	//Pointers to a vector of dimendion N that represent forces acting on lattice points, in initialise()	      (SISTEMA)
	delete []forcex; delete []forcey; delete []forcez;
	//Pointers to a vector of dimension N that is the neighbours of each lattice point, in initialise()
	delete []dd1; delete []dd2; delete []dd3; 
	delete []dd4; delete []dd5; delete []dd6; 
	delete []dda; delete []ddb; delete []ddc; 
	delete []ddd; delete []dde; delete []ddf;
	delete []ddg; delete []ddh; delete []ddi;
	delete []ddj; delete []ddk; delete []ddl;
  	//Pointers to a vector of dimension N that is lattice ID for topolocially patterned substrate, defined in initialise()
	delete []mask;

	delete []maskGlobal;
	delete []pGlobal;

	
	cout << "Process " << rank << ": Destructor: variables space disallocated" << endl;
}

