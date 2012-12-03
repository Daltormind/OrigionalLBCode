//writeYPlanXVelocityFile(Ylayer).cpp: for a certain Y-layer it writes the x and z velocity components into matlab m-file: uxt(:,:)= [...]
//filename: uxYytT.m where T = number of iterations and y is the Y-layer number
//Lisa, 28 July 2011
// Changed to write velocities in each Z layer instead, Mat 19 Oct 2012

#include "wet.h"



void wet::writeZPlanXVelocityFile(const int iz)
{
	
	int j=iz-1;
	int l,i,h;
	 
	
	//cout << " Process " << rank << ": check1 " << endl;
	MPI_Gather(&(uxs[k1]),k2-k1,MPI_DOUBLE,uGlobal,k2-k1,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);
	//cout << " Process " << rank << ": check2 " << endl;
	//cout << "k " << k <<  " uGlobal " << uGlobal << endl; 

	if (rank == ROOT) {
		
		char filenamec[20];
		string filename;
        
        snprintf(filenamec,20, "/duxZ%dt%ld.m" ,iz,t);
		filename=folder+filenamec;
        
        ofstream file(filename.c_str());
		
		file.precision(4);
		file << "uxZ"<< iz << "t" << t << "=[" << endl;
		for(h = 0 ; h < LY ; h++) {
			for(i = 0 ; i < LX ; i++) {   
				l= h + j*LZ + i*LY*LZ;	
				file << uGlobal[l] << " ";
				//cout << "at t = " << t << " uxGl at x= " << i << " y= " << j << " z= " << h << " is " << uGlobal[l] << endl ;
			}
			file << endl;
		}
		file << "];" << endl;
		file.close();
		
		/*for (k=0; k<N; k++) {
			
			xk=int(k/(float) (LZ*LY));
			yk=int((k-xk*LZ*LY)/(float) LZ);
			zk=k-xk*LZ*LY-yk*LZ;
			
			cout << "AT T = " << t << " uxGl at x= " << xk << " y= " << yk << " z= " << zk << " is " << uGlobal[k] << endl ;
		}*/

	}
	/*
	if (rank == ROOT ) {
		
		char filename[30];
		sprintf(filename, "./densitydata/dvelx%ld.m", t);	
		ofstream file(filename);
		
		file.precision(4);
		for( h = 0 ; h < LZ ; h++) 
		{   
			file << "velx" << t << "(:,:," << h+1 << ")=[" << endl;
			for( i = 0 ; i < LX ; i++) 
			{
				for( j = 0 ; j < LY ; j++) 
				{
					k = h + j*LZ + i*LY*LZ;
					if( maskGlobal[k]==28) 
						file << 0.0 << " ";  
					else
						file << uGlobal[k] << " ";
					
				}
				file << endl;
			}
			file <<"];" << endl;
		}
		
		file.close();
		
	}
	*/
	
	
	
	
	
	MPI_Gather(&(uys[k1]),k2-k1,MPI_DOUBLE,uGlobal,k2-k1,MPI_DOUBLE,0,MPI_COMM_WORLD);

	if (rank == ROOT) {
		
		char filename[20];
		snprintf(filename,20, "/duyZ%dt%ld.m", iz,t);
		string filename1;
        filename1=folder+filename;
        ofstream file(filename1.c_str());
		file.precision(4);
		file << "uyZ"<< iz << "t" << t << "=[" << endl;

		for(h = 0 ; h < LY ; h++) {
			for(i = 0 ; i < LX ; i++) {   
				l= h + j*LZ + i*LY*LZ;	
				file << uGlobal[l] << " ";
				//cout << "at t = " << t << " uxGl at x= " << i << " y= " << j << " z= " << h << " is " << uGlobal[l] << endl ;
			}
			file << endl;
		}
		file << "];" << endl;
		file.close();
		
	}
	
	// write layer in y direction at x=0:
	/*
	int ix = LX/2;
	int m=ix-1;
	
	MPI_Gather(&(uys[k1]),k2-k1,MPI_DOUBLE,uGlobal,k2-k1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	
	if (rank == ROOT) {
		
		char filename[30];
		sprintf(filename, "./densitydata/duyX%dt%ld.m",ix,t);
	
		ofstream file(filename);
		
		file.precision(4);
		file << "uyX"<< ix << "t" << t << "=[" << endl;
		
		for(h = 0 ; h < LZ ; h++) {
			for(i = 0 ; i < LY ; i++) {   
				l= h + i*LZ + m*LY*LZ;	
				file << uGlobal[l] << " ";
				//cout << "at t = " << t << " uxGl at x= " << i << " y= " << j << " z= " << h << " is " << uGlobal[l] << endl ;
			}
			file << endl;
		}
		file << "];" << endl;
		file.close();
		
	}
	
	MPI_Gather(&(uzs[k1]),k2-k1,MPI_DOUBLE,uGlobal,k2-k1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	
	if (rank == ROOT) {
		
		char filename[30];
		sprintf(filename, "./densitydata/duzX%dt%ld.m",ix,t);
		
		ofstream file(filename);
		
		file.precision(4);
		file << "uzX"<< ix << "t" << t << "=[" << endl;
		
		for(h = 0 ; h < LZ ; h++) {
			for(i = 0 ; i < LY ; i++) {   
				l= h + i*LZ + m*LY*LZ;	
				file << uGlobal[l] << " ";
				//cout << "at t = " << t << " uxGl at x= " << i << " y= " << j << " z= " << h << " is " << uGlobal[l] << endl ;
			}
			file << endl;
		}
		file << "];" << endl;
		file.close();
		
	}*/
	
	
	
	/*
	if (rank == ROOT) {
		
		char filename[30];
		sprintf(filename, "./densitydata/dvelz%ld.m", t);	
		ofstream file(filename);
		
		file.precision(4);
		for( h = 0 ; h < LZ ; h++) 
		{   
			file << "velz" << t << "(:,:," << h+1 << ")=[" << endl;
			for( i = 0 ; i < LX ; i++) 
			{
				for( j = 0 ; j < LY ; j++) 
				{
					k = h + j*LZ + i*LY*LZ;
					if( maskGlobal[k]==28) 
						file << 0.0 << " ";  
					else
						file << uGlobal[k] << " ";
					
				}
				file << endl;
			}
			file <<"];" << endl;
		}
		
		
		file.close();
		
	}
	
	MPI_Gather(&(uys[k1]),k2-k1,MPI_DOUBLE,uGlobal,k2-k1,MPI_DOUBLE,0,MPI_COMM_WORLD);

	if (rank == ROOT) {
		
				
		char filename[30];
		sprintf(filename, "./densitydata/dvely%ld.m", t);	
		ofstream file(filename);
		
		file.precision(4);
		for( h = 0 ; h < LZ ; h++) 
		{   
			file << "vely" << t << "(:,:," << h+1 << ")=[" << endl;
			for( i = 0 ; i < LX ; i++) 
			{
				for( j = 0 ; j < LY ; j++) 
				{
					k = h + j*LZ + i*LY*LZ;
					if( maskGlobal[k]==28) 
						file << 0.0 << " ";  
					else
						file << uGlobal[k] << " ";
					
				}
				file << endl;
			}
			file <<"];" << endl;
		}
		
		
		file.close();
		
	}
	*/
}

/*
 String fileName("uxY");
 fileName.concat((int) j);
 fileName.concat("t");
 fileName.concat((int) t);
 fileName.concat(".m");
 ofstream file(fileName.get());
 */

/*
 String fileName("uzY");
 fileName.concat((int) j);
 fileName.concat("t");
 fileName.concat((int) t);
 fileName.concat(".m");
 ofstream file(fileName.get());
 */