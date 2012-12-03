// writehandlefile.cpp write a file which can be run as a matlab script to give parameters for data analysis. 

#include "wet.h"


void wet::writehandlefile()
{

	char filename[5];
	string filenames;
	snprintf(filename,5,"/p.m");
	filenames=folder+filename;
	ofstream file(filenames.c_str());
	

	file << "writestep=" << writeStep << ";" << endl;
	file << "nbEqStep=" << nbEqStep << ";" << endl;
	file << "infostep=" << infoStep << ";" << endl;
	file << "folder=" << "'" << folder << "';" << endl;
	

}
	