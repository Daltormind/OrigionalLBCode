#include "wet.h"

void wet::makematrix(void)
{
	double tau;
	double S[19];
	double Z[19][19];
	int mag[19];
	int i, j, m;
	
	cout << "Process " << rank << ": creation \"evolution matrix\"..." << endl;

	for(j = 0; j < 19; j++)
	{ 
		mag[j] = 0; 
		for(i = 0; i < 19; i++) 
			mag[j] += M[j][i]*M[j][i];
	}
  
	for(l = 0; l < NRATIOFLUID; l++)
	{
		tau = taugas + double(l)*(tauliquid-taugas)/(NRATIOFLUID-1);
		S[0] = 0.0; S[1] = 1.0; S[2] = 1.0; S[3] = 0.0; S[4] = 1.0; S[5] = 0.0; S[6] = 1.0; 
		S[7] = 0.0; S[8] = 1.0; S[9] = 1.0/tau; S[10] = 1.0; S[11] = 1.0/tau; S[12] = 1.0; 
		S[13] = 1.0/tau; S[14] = 1.0/tau; S[15] = 1.0/tau; S[16] = 1.0; S[17] = 1.0; S[18] = 1.0;
		
		for(i=0; i<19; i++)
		{
			for(j=0; j<19; j++)
			{
				Z[i][j]=0;
				for(m=0; m<19; m++)
					Z[i][j]+= M[m][i]*S[m]*M[m][j]/mag[m];
					
			}
		}
		 
    		  
		M00[l]=Z[0][0];  M01[l] = Z[0][1];  M02[l] = Z[0][2];  M03[l] = Z[0][3];  M04[l] = Z[0][4];  M05[l] = Z[0][5];  M06[l] = Z[0][6]; 
		M10[l]=Z[1][0];  M11[l] = Z[1][1];  M12[l] = Z[1][2];  M13[l] = Z[1][3];  M14[l] = Z[1][4];  M15[l] = Z[1][5];  M16[l] = Z[1][6]; 
		M20[l] = Z[2][0];  M21[l] = Z[2][1];  M22[l] = Z[2][2];  M23[l] = Z[2][3];  M24[l] = Z[2][4];  M25[l] = Z[2][5];  M26[l] = Z[2][6]; 
		M30[l] = Z[3][0];  M31[l] = Z[3][1];  M32[l] = Z[3][2];  M33[l] = Z[3][3];  M34[l] = Z[3][4];  M35[l] = Z[3][5];  M36[l] = Z[3][6]; 
		M40[l] = Z[4][0];  M41[l] = Z[4][1];  M42[l] = Z[4][2];  M43[l] = Z[4][3];  M44[l] = Z[4][4];  M45[l] = Z[4][5];  M46[l] = Z[4][6]; 
		M50[l] = Z[5][0];  M51[l] = Z[5][1];  M52[l] = Z[5][2];  M53[l] = Z[5][3];  M54[l] = Z[5][4];  M55[l] = Z[5][5];  M56[l] = Z[5][6]; 
		M60[l] = Z[6][0];  M61[l] = Z[6][1];  M62[l] = Z[6][2];  M63[l] = Z[6][3];  M64[l] = Z[6][4];  M65[l] = Z[6][5];  M66[l] = Z[6][6]; 
		Ma0[l] = Z[7][0];  Ma1[l] = Z[7][1];  Ma2[l] = Z[7][2];  Ma3[l] = Z[7][3];  Ma4[l] = Z[7][4];  Ma5[l] = Z[7][5];  Ma6[l] = Z[7][6]; 
		Mb0[l] = Z[8][0];  Mb1[l] = Z[8][1];  Mb2[l] = Z[8][2];  Mb3[l] = Z[8][3];  Mb4[l] = Z[8][4];  Mb5[l] = Z[8][5];  Mb6[l] = Z[8][6]; 
		Mc0[l] = Z[9][0];  Mc1[l] = Z[9][1];  Mc2[l] = Z[9][2];  Mc3[l] = Z[9][3];  Mc4[l] = Z[9][4];  Mc5[l] = Z[9][5];  Mc6[l] = Z[9][6]; 
		Md0[l] = Z[10][0]; Md1[l] = Z[10][1]; Md2[l] = Z[10][2]; Md3[l] = Z[10][3]; Md4[l] = Z[10][4]; Md5[l] = Z[10][5]; Md6[l] = Z[10][6]; 
		Me0[l] = Z[11][0]; Me1[l] = Z[11][1]; Me2[l] = Z[11][2]; Me3[l] = Z[11][3]; Me4[l] = Z[11][4]; Me5[l] = Z[11][5]; Me6[l] = Z[11][6]; 
		Mf0[l] = Z[12][0]; Mf1[l] = Z[12][1]; Mf2[l] = Z[12][2]; Mf3[l] = Z[12][3]; Mf4[l] = Z[12][4]; Mf5[l] = Z[12][5]; Mf6[l] = Z[12][6]; 
		Mg0[l] = Z[13][0]; Mg1[l] = Z[13][1]; Mg2[l] = Z[13][2]; Mg3[l] = Z[13][3]; Mg4[l] = Z[13][4]; Mg5[l] = Z[13][5]; Mg6[l] = Z[13][6]; 
		Mh0[l] = Z[14][0]; Mh1[l] = Z[14][1]; Mh2[l] = Z[14][2]; Mh3[l] = Z[14][3]; Mh4[l] = Z[14][4]; Mh5[l] = Z[14][5]; Mh6[l] = Z[14][6]; 
		Mi0[l] = Z[15][0]; Mi1[l] = Z[15][1]; Mi2[l] = Z[15][2]; Mi3[l] = Z[15][3]; Mi4[l] = Z[15][4]; Mi5[l] = Z[15][5]; Mi6[l] = Z[15][6]; 
		Mj0[l] = Z[16][0]; Mj1[l] = Z[16][1]; Mj2[l] = Z[16][2]; Mj3[l] = Z[16][3]; Mj4[l] = Z[16][4]; Mj5[l] = Z[16][5]; Mj6[l] = Z[16][6]; 
		Mk0[l] = Z[17][0]; Mk1[l] = Z[17][1]; Mk2[l] = Z[17][2]; Mk3[l] = Z[17][3]; Mk4[l] = Z[17][4]; Mk5[l] = Z[17][5]; Mk6[l] = Z[17][6]; 
		Ml0[l] = Z[18][0]; Ml1[l] = Z[18][1]; Ml2[l] = Z[18][2]; Ml3[l] = Z[18][3]; Ml4[l] = Z[18][4]; Ml5[l] = Z[18][5]; Ml6[l] = Z[18][6]; 
    
		M0a[l] = Z[0][7];  M0b[l] = Z[0][8];  M0c[l] = Z[0][9];  M0d[l] = Z[0][10];  M0e[l] = Z[0][11];  M0f[l] = Z[0][12];
		M1a[l] = Z[1][7];  M1b[l] = Z[1][8];  M1c[l] = Z[1][9];  M1d[l] = Z[1][10];  M1e[l] = Z[1][11];  M1f[l] = Z[1][12]; 
		M2a[l] = Z[2][7];  M2b[l] = Z[2][8];  M2c[l] = Z[2][9];  M2d[l] = Z[2][10];  M2e[l] = Z[2][11];  M2f[l] = Z[2][12]; 
		M3a[l] = Z[3][7];  M3b[l] = Z[3][8];  M3c[l] = Z[3][9];  M3d[l] = Z[3][10];  M3e[l] = Z[3][11];  M3f[l] = Z[3][12];
		M4a[l] = Z[4][7];  M4b[l] = Z[4][8];  M4c[l] = Z[4][9];  M4d[l] = Z[4][10];  M4e[l] = Z[4][11];  M4f[l] = Z[4][12]; 
		M5a[l] = Z[5][7];  M5b[l] = Z[5][8];  M5c[l] = Z[5][9];  M5d[l] = Z[5][10];  M5e[l] = Z[5][11];  M5f[l] = Z[5][12]; 
		M6a[l] = Z[6][7];  M6b[l] = Z[6][8];  M6c[l] = Z[6][9];  M6d[l] = Z[6][10];  M6e[l] = Z[6][11];  M6f[l] = Z[6][12]; 
		Maa[l] = Z[7][7];  Mab[l] = Z[7][8];  Mac[l] = Z[7][9];  Mad[l] = Z[7][10];  Mae[l] = Z[7][11];  Maf[l] = Z[7][12];
		Mba[l] = Z[8][7];  Mbb[l] = Z[8][8];  Mbc[l] = Z[8][9];  Mbd[l] = Z[8][10];  Mbe[l] = Z[8][11];  Mbf[l] = Z[8][12]; 
		Mca[l] = Z[9][7];  Mcb[l] = Z[9][8];  Mcc[l] = Z[9][9];  Mcd[l] = Z[9][10];  Mce[l] = Z[9][11];  Mcf[l] = Z[9][12];
		Mda[l] = Z[10][7]; Mdb[l] = Z[10][8]; Mdc[l] = Z[10][9]; Mdd[l] = Z[10][10]; Mde[l] = Z[10][11]; Mdf[l] = Z[10][12]; 
		Mea[l] = Z[11][7]; Meb[l] = Z[11][8]; Mec[l] = Z[11][9]; Med[l] = Z[11][10]; Mee[l] = Z[11][11]; Mef[l] = Z[11][12];
		Mfa[l] = Z[12][7]; Mfb[l] = Z[12][8]; Mfc[l] = Z[12][9]; Mfd[l] = Z[12][10]; Mfe[l] = Z[12][11]; Mff[l] = Z[12][12]; 
		Mga[l] = Z[13][7]; Mgb[l] = Z[13][8]; Mgc[l] = Z[13][9]; Mgd[l] = Z[13][10]; Mge[l] = Z[13][11]; Mgf[l] = Z[13][12]; 
		Mha[l] = Z[14][7]; Mhb[l] = Z[14][8]; Mhc[l] = Z[14][9]; Mhd[l] = Z[14][10]; Mhe[l] = Z[14][11]; Mhf[l] = Z[14][12]; 
		Mia[l] = Z[15][7]; Mib[l] = Z[15][8]; Mic[l] = Z[15][9]; Mid[l] = Z[15][10]; Mie[l] = Z[15][11]; Mif[l] = Z[15][12]; 
		Mja[l] = Z[16][7]; Mjb[l] = Z[16][8]; Mjc[l] = Z[16][9]; Mjd[l] = Z[16][10]; Mje[l] = Z[16][11]; Mjf[l] = Z[16][12]; 
		Mka[l] = Z[17][7]; Mkb[l] = Z[17][8]; Mkc[l] = Z[17][9]; Mkd[l] = Z[17][10]; Mke[l] = Z[17][11]; Mkf[l] = Z[17][12]; 
		Mla[l] = Z[18][7]; Mlb[l] = Z[18][8]; Mlc[l] = Z[18][9]; Mld[l] = Z[18][10]; Mle[l] = Z[18][11]; Mlf[l] = Z[18][12];
    
		M0g[l] = Z[0][13];  M0h[l] = Z[0][14];  M0i[l] = Z[0][15];  M0j[l] = Z[0][16];  M0k[l] = Z[0][17];  M0l[l] = Z[0][18];
		M1g[l] = Z[1][13];  M1h[l] = Z[1][14];  M1i[l] = Z[1][15];  M1j[l] = Z[1][16];  M1k[l] = Z[1][17];  M1l[l] = Z[1][18]; 
		M2g[l] = Z[2][13];  M2h[l] = Z[2][14];  M2i[l] = Z[2][15];  M2j[l] = Z[2][16];  M2k[l] = Z[2][17];  M2l[l] = Z[2][18]; 
		M3g[l] = Z[3][13];  M3h[l] = Z[3][14];  M3i[l] = Z[3][15];  M3j[l] = Z[3][16];  M3k[l] = Z[3][17];  M3l[l] = Z[3][18];
		M4g[l] = Z[4][13];  M4h[l] = Z[4][14];  M4i[l] = Z[4][15];  M4j[l] = Z[4][16];  M4k[l] = Z[4][17];  M4l[l] = Z[4][18]; 
		M5g[l] = Z[5][13];  M5h[l] = Z[5][14];  M5i[l] = Z[5][15];  M5j[l] = Z[5][16];  M5k[l] = Z[5][17];  M5l[l] = Z[5][18]; 
		M6g[l] = Z[6][13];  M6h[l] = Z[6][14];  M6i[l] = Z[6][15];  M6j[l] = Z[6][16];  M6k[l] = Z[6][17];  M6l[l] = Z[6][18]; 
		Mag[l] = Z[7][13];  Mah[l] = Z[7][14];  Mai[l] = Z[7][15];  Maj[l] = Z[7][16];  Mak[l] = Z[7][17];  Mal[l] = Z[7][18];
		Mbg[l] = Z[8][13];  Mbh[l] = Z[8][14];  Mbi[l] = Z[8][15];  Mbj[l] = Z[8][16];  Mbk[l] = Z[8][17];  Mbl[l] = Z[8][18]; 
		Mcg[l] = Z[9][13];  Mch[l] = Z[9][14];  Mci[l] = Z[9][15];  Mcj[l] = Z[9][16];  Mck[l] = Z[9][17];  Mcl[l] = Z[9][18];
		Mdg[l] = Z[10][13]; Mdh[l] = Z[10][14]; Mdi[l] = Z[10][15]; Mdj[l] = Z[10][16]; Mdk[l] = Z[10][17]; Mdl[l] = Z[10][18]; 
		Meg[l] = Z[11][13]; Meh[l] = Z[11][14]; Mei[l] = Z[11][15]; Mej[l] = Z[11][16]; Mek[l] = Z[11][17]; Mel[l] = Z[11][18];
		Mfg[l] = Z[12][13]; Mfh[l] = Z[12][14]; Mfi[l] = Z[12][15]; Mfj[l] = Z[12][16]; Mfk[l] = Z[12][17]; Mfl[l] = Z[12][18]; 
		Mgg[l] = Z[13][13]; Mgh[l] = Z[13][14]; Mgi[l] = Z[13][15]; Mgj[l] = Z[13][16]; Mgk[l] = Z[13][17]; Mgl[l] = Z[13][18]; 
		Mhg[l] = Z[14][13]; Mhh[l] = Z[14][14]; Mhi[l] = Z[14][15]; Mhj[l] = Z[14][16]; Mhk[l] = Z[14][17]; Mhl[l] = Z[14][18]; 
		Mig[l] = Z[15][13]; Mih[l] = Z[15][14]; Mii[l] = Z[15][15]; Mij[l] = Z[15][16]; Mik[l] = Z[15][17]; Mil[l] = Z[15][18]; 
		Mjg[l] = Z[16][13]; Mjh[l] = Z[16][14]; Mji[l] = Z[16][15]; Mjj[l] = Z[16][16]; Mjk[l] = Z[16][17]; Mjl[l] = Z[16][18]; 
		Mkg[l] = Z[17][13]; Mkh[l] = Z[17][14]; Mki[l] = Z[17][15]; Mkj[l] = Z[17][16]; Mkk[l] = Z[17][17]; Mkl[l] = Z[17][18]; 
		Mlg[l] = Z[18][13]; Mlh[l] = Z[18][14]; Mli[l] = Z[18][15]; Mlj[l] = Z[18][16]; Mlk[l] = Z[18][17]; Mll[l] = Z[18][18];
	}
	
	cout << "Process " << rank << ": matrix created." << endl;
}
