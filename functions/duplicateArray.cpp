//duplicateArray.cpp: takes an array of size LXc/2*LY*LZ = Nc/2 and copies it over to a new array of N = Nc = LXc*LY*LZ
//such that the new part is the mirror image of the existing part
//Lisa, 02 Aug 2011

#include "wet.h"

double* wet::duplicateArray(double *oldArr)
{	
	int kdup;
	double *newArr;
	
	newArr = new double[processN]; 
	
	//MPI_Barrier(MPI_COMM_WORLD);
		
	if (rank < size/2){          //first half of the box, drop just copied over
		
		
		for (k=k1 ; k< k2; k++) {
			
			
			kdup = k-k1 + rank * (k2-k1);
			
			/*
			xkdup =int(kdup/(float) (LZ*LY));
			ykdup =int((kdup-xkdup*LZ*LY)/(float) LZ);
			zkdup =kdup-xkdup*LZ*LY-ykdup*LZ;
			*/
			
			newArr[k] = oldArr[kdup];
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
			
			/*
			xkdup =int(kdup/(float) (LZ*LY));
			ykdup =int((kdup-xkdup*LZ*LY)/(float) LZ);
			zkdup =kdup-xkdup*LZ*LY-ykdup*LZ;
			xkdup =xkdup+rank*LX/size-1%LX;
			 */
			
			newArr[kdup] = oldArr[k];
			
			//computeCoordinate(k);
			//cout << "<< Process " << rank <<" k= "<< kdup << " xk yk zk " << xkdup << " " << ykdup << " " << zkdup <<" " << "p[k] = " << p[kdup] <<"; kdup= "<< k << " xk yk zk " << xk << " " << yk << " " << zk << " " << "pG[kdup] = " << pGlobal[k] << endl;
			
			
		}
		
	}
	
	return newArr;	
}
	
/*double* wet::duplicateArray(int *oldArr)
{	
	int kdup;
	double *newArr;
	
	newArr = new double[N];
		
	if (rank < size/2){          //first half of the box, drop just copied over
		
		
		for (k=k1 ; k< k2; k++) {
			
			
			kdup = k-k1 + rank * (k2-k1);
			
			newArr[k] = double(oldArr[kdup]);
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
			
			newArr[kdup] = double(oldArr[k]);
			
			//computeCoordinate(k);
			//cout << "<< Process " << rank <<" k= "<< kdup << " xk yk zk " << xkdup << " " << ykdup << " " << zkdup <<" " << "p[k] = " << p[kdup] <<"; kdup= "<< k << " xk yk zk " << xk << " " << yk << " " << zk << " " << "pG[kdup] = " << pGlobal[k] << endl;
			
			
		}
		
	}
	
	return newArr;	
}
*/
