//Parallel process debug function

	char filenamefn1[30];
		sprintf(filenamefn1, "./data/process%d.fn1.t%ld.dat",rank, t);	
	char filenamegn1[30];
		sprintf(filenamegn1, "./data/process%d.gn1.t%ld.dat",rank, t);	
	char filenamefn2[30];
		sprintf(filenamefn2, "./data/process%d.fn2.t%ld.dat",rank, t);	
	char filenamegn2[30];
		sprintf(filenamegn2, "./data/process%d.gn2.t%ld.dat",rank, t);	
	char filenamefna[30];
		sprintf(filenamefna, "./data/process%d.fna.t%ld.dat",rank, t);	
	char filenamegna[30];
		sprintf(filenamegna, "./data/process%d.gna.t%ld.dat",rank, t);	
	char filenamefnd[30];
		sprintf(filenamefnd, "./data/process%d.fnd.t%ld.dat",rank, t);	
	char filenamegnd[30];
		sprintf(filenamegnd, "./data/process%d.gnd.t%ld.dat",rank, t);	
	char filenamefnb[30];
		sprintf(filenamefnb, "./data/process%d.fnb.t%ld.dat",rank, t);	
	char filenamegnb[30];
		sprintf(filenamegnb, "./data/process%d.gnb.t%ld.dat",rank, t);	
	char filenamefnc[30];
		sprintf(filenamefnc, "./data/process%d.fnc.t%ld.dat",rank, t);	
	char filenamegnc[30];
		sprintf(filenamegnc, "./data/process%d.gnc.t%ld.dat",rank, t);	
	char filenamefni[30];
		sprintf(filenamefni, "./data/process%d.fni.t%ld.dat",rank, t);	
	char filenamegni[30];
		sprintf(filenamegni, "./data/process%d.gni.t%ld.dat",rank, t);
	char filenamefnl[30];
		sprintf(filenamefnl, "./data/process%d.fnl.t%ld.dat",rank, t);	
	char filenamegnl[30];
		sprintf(filenamegnl, "./data/process%d.gnl.t%ld.dat",rank, t);	
	char filenamefnj[30];
		sprintf(filenamefnj, "./data/process%d.fnj.t%ld.dat",rank, t);	
	char filenamegnj[30];
		sprintf(filenamegnj, "./data/process%d.gnj.t%ld.dat",rank, t);
	char filenamefnk[30];
		sprintf(filenamefnk, "./data/process%d.fnk.t%ld.dat",rank, t);	
	char filenamegnk[30];
		sprintf(filenamegnk, "./data/process%d.gnk.t%ld.dat",rank, t);	
	

	

		ofstream filefn1(filenamefn1);
		for(k=0; k<processN; k++)
		{
			filefn1 << fn1[k] << "	";
			if((k+1)%LZ==0 && k!=0)
			filefn1 << endl;
			if((k+1)%k1==0 && k!=0)
			filefn1 << endl << endl << endl;
		}
		filefn1.close();

		ofstream filegn1(filenamegn1);
		for(k=0; k<processN; k++)
		{
			filegn1 << gn1[k] << "	";
			if((k+1)%LZ==0)
			filegn1 << endl;
			if((k+1)%k1==0)
			filegn1 << endl << endl << endl;
		}
		filegn1.close();
	
		ofstream filefn2(filenamefn2);
		for(k=0; k<processN; k++)
		{
			filefn2 << fn2[k] << "	";
			if((k+1)%LZ==0 && k!=0)
			filefn2 << endl;
			if((k+1)%k1==0 && k!=0)
			filefn2 << endl << endl << endl;
		}
		filefn2.close();

		ofstream filegn2(filenamegn2);
		for(k=0; k<processN; k++)
		{
			filegn2 << gn2[k] << "	";
			if((k+1)%LZ==0)
			filegn2 << endl;
			if((k+1)%k1==0)
			filegn2 << endl << endl << endl;
		}
		filegn2.close();

		ofstream filefna(filenamefna);
		for(k=0; k<processN; k++)
		{
			filefna << fna[k] << "	";
			if((k+1)%LZ==0 && k!=0)
			filefna << endl;
			if((k+1)%k1==0 && k!=0)
			filefna << endl << endl << endl;
		}
		filefna.close();

		ofstream filegna(filenamegna);
		for(k=0; k<processN; k++)
		{
			filegna << gna[k] << "	";
			if((k+1)%LZ==0)
			filegna << endl;
			if((k+1)%k1==0)
			filegna << endl << endl << endl;
		}
		filegna.close();		

		ofstream filefnd(filenamefnd);
		for(k=0; k<processN; k++)
		{
			filefnd << fnd[k] << "	";
			if((k+1)%LZ==0 && k!=0)
			filefnd << endl;
			if((k+1)%k1==0 && k!=0)
			filefnd << endl << endl << endl;
		}
		filefnd.close();

		ofstream filegnd(filenamegnd);
		for(k=0; k<processN; k++)
		{
			filegnd << gnd[k] << "	";
			if((k+1)%LZ==0)
			filegnd << endl;
			if((k+1)%k1==0)
			filegnd << endl << endl << endl;
		}
		filegnd.close();

				ofstream filefnb(filenamefnb);
		for(k=0; k<processN; k++)
		{
			filefnb << fnb[k] << "	";
			if((k+1)%LZ==0 && k!=0)
			filefnb << endl;
			if((k+1)%k1==0 && k!=0)
			filefnb << endl << endl << endl;
		}
		filefnb.close();

		ofstream filegnb(filenamegnb);
		for(k=0; k<processN; k++)
		{
			filegnb << gnb[k] << "	";
			if((k+1)%LZ==0)
			filegnb << endl;
			if((k+1)%k1==0)
			filegnb << endl << endl << endl;
		}
		filegnb.close();

		ofstream filefnc(filenamefnc);
		for(k=0; k<processN; k++)
		{
			filefnc << fnc[k] << "	";
			if((k+1)%LZ==0 && k!=0)
			filefnc << endl;
			if((k+1)%k1==0 && k!=0)
			filefnc << endl << endl << endl;
		}
		filefnc.close();

		ofstream filegnc(filenamegnc);
		for(k=0; k<processN; k++)
		{
			filegnc << gnc[k] << "	";
			if((k+1)%LZ==0)
			filegnc << endl;
			if((k+1)%k1==0)
			filegnc << endl << endl << endl;
		}
		filegnc.close();

	
		ofstream filefni(filenamefni);
		for(k=0; k<processN; k++)
		{
			filefni << fni[k] << "	";
			if((k+1)%LZ==0 && k!=0)
			filefni << endl;
			if((k+1)%k1==0 && k!=0)
			filefni << endl << endl << endl;
		}
		filefni.close();

		ofstream filegni(filenamegni);
		for(k=0; k<processN; k++)
		{
			filegni << gni[k] << "	";
			if((k+1)%LZ==0)
			filegni << endl;
			if((k+1)%k1==0)
			filegni << endl << endl << endl;
		}
		filegni.close();

		ofstream filefnl(filenamefnl);
		for(k=0; k<processN; k++)
		{
			filefnl << fnl[k] << "	";
			if((k+1)%LZ==0 && k!=0)
			filefnl << endl;
			if((k+1)%k1==0 && k!=0)
			filefnl << endl << endl << endl;
		}
		filefnl.close();

		ofstream filegnl(filenamegnl);
		for(k=0; k<processN; k++)
		{
			filegnl << gnl[k] << "	";
			if((k+1)%LZ==0)
			filegnl << endl;
			if((k+1)%k1==0)
			filegnl << endl << endl << endl;
		}
		filegnl.close();

		ofstream filefnj(filenamefnj);
		for(k=0; k<processN; k++)
		{
			filefnj << fnj[k] << "	";
			if((k+1)%LZ==0 && k!=0)
			filefnj << endl;
			if((k+1)%k1==0 && k!=0)
			filefnj << endl << endl << endl;
		}
		filefnj.close();

		ofstream filegnj(filenamegnj);
		for(k=0; k<processN; k++)
		{
			filegnj << gnj[k] << "	";
			if((k+1)%LZ==0)
			filegnj << endl;
			if((k+1)%k1==0)
			filegnj << endl << endl << endl;
		}
		filegnj.close();

		ofstream filefnk(filenamefnk);
		for(k=0; k<processN; k++)
		{
			filefnk << fnk[k] << "	";
			if((k+1)%LZ==0 && k!=0)
			filefnk << endl;
			if((k+1)%k1==0 && k!=0)
			filefnk << endl << endl << endl;
		}
		filefnk.close();

		ofstream filegnk(filenamegnk);
		for(k=0; k<processN; k++)
		{
			filegnk << gnk[k] << "	";
			if((k+1)%LZ==0)
			filegnk << endl;
			if((k+1)%k1==0)
			filegnk << endl << endl << endl;
		}
		filegnk.close();

//write file .dat for each process
	char filename1[30];
		sprintf(filename1, "./data/process%d.n%ld.dat",rank, t);	
	char filename2[30];
		sprintf(filename2, "./data/process%d.p%ld.dat",rank, t);	
	
		ofstream filea(filename1);
		for(k=0; k<processN; k++)
		{
			filea << n[k] << "	";
			if((k+1)%LZ==0 && k!=0)
			filea << endl;
			if((k+1)%k1==0 && k!=0)
			filea << endl << endl << endl;
		}
		filea.close();

		ofstream fileb(filename2);
		for(k=0; k<processN; k++)
		{
			fileb << p[k] << "	";
			if((k+1)%LZ==0)
			fileb << endl;
			if((k+1)%k1==0)
			fileb << endl << endl << endl;
		}
		fileb.close();

//write phi derivatives for each process
	cout << "Process " << rank << ": dp_dx(" << k <<")= " << dp_dx << endl;
	cout << "Process " << rank << ": dp_dy(" << k <<")= " << dp_dy << endl;
	cout << "Process " << rank << ": dp_dz(" << k <<")= " << dp_dz << endl;

