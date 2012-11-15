//collision.cpp: function that produces collision between nodes
//Marco, 30 June 2011

#include "wet.h"

void wet::collision(void)
{ 

	//if(t%infoStep==0)
	//	cout << "Process " << rank << ": colliding..." << endl;
	
	for( k = k1 ; k < k2 ; k++)
	{
		if(mask[k]!=28)
		{

			equilibrium();
			l = int((p[k]+1)*double(NRATIOFLUID)/2.0);
			if(l < 0) 
				l = 0; 
			if(l > NRATIOFLUID-1) 
				l = NRATIOFLUID-1;
   
			df0 = fe0-ff0[k]; df1 = fe1-ff1[k]; df2 = fe2-ff2[k]; df3 = fe3-ff3[k]; df4 = fe4-ff4[k]; df5 = fe5-ff5[k]; df6 = fe6-ff6[k]; 
			dfa = fea-ffa[k]; dfb = feb-ffb[k]; dfc = fec-ffc[k]; dfd = fed-ffd[k]; dfe = fee-ffe[k]; dff = fef-fff[k]; dfg = feg-ffg[k]; dfh = feh-ffh[k]; dfi = fei-ffi[k]; dfj = fej-ffj[k]; dfk = fek-ffk[k]; dfl = fel-ffl[k];

			fn0[k]  = ff0[k] + df0*M00[l] + df1*M01[l] + df2*M02[l] + df3*M03[l] + df4*M04[l] + df5*M05[l] + df6*M06[l] + dfa*M0a[l] + dfb*M0b[l] + dfc*M0c[l] + dfd*M0d[l] + dfe*M0e[l] + dff*M0f[l] + dfg*M0g[l] + dfh*M0h[l] + dfi*M0i[l] + dfj*M0j[l] + dfk*M0k[l] + dfl*M0l[l] + force0; 
			fn1[k]  = ff1[k] + df0*M10[l] + df1*M11[l] + df2*M12[l] + df3*M13[l] + df4*M14[l] + df5*M15[l] + df6*M16[l] + dfa*M1a[l] + dfb*M1b[l] + dfc*M1c[l] + dfd*M1d[l] + dfe*M1e[l] + dff*M1f[l] + dfg*M1g[l] + dfh*M1h[l] + dfi*M1i[l] + dfj*M1j[l] + dfk*M1k[l] + dfl*M1l[l] + force1; 
			fn2[k]  = ff2[k] + df0*M20[l] + df1*M21[l] + df2*M22[l] + df3*M23[l] + df4*M24[l] + df5*M25[l] + df6*M26[l] + dfa*M2a[l] + dfb*M2b[l] + dfc*M2c[l] + dfd*M2d[l] + dfe*M2e[l] + dff*M2f[l] + dfg*M2g[l] + dfh*M2h[l] + dfi*M2i[l] + dfj*M2j[l] + dfk*M2k[l] + dfl*M2l[l] + force2; 
			fn3[k]  = ff3[k] + df0*M30[l] + df1*M31[l] + df2*M32[l] + df3*M33[l] + df4*M34[l] + df5*M35[l] + df6*M36[l] + dfa*M3a[l] + dfb*M3b[l] + dfc*M3c[l] + dfd*M3d[l] + dfe*M3e[l] + dff*M3f[l] + dfg*M3g[l] + dfh*M3h[l] + dfi*M3i[l] + dfj*M3j[l] + dfk*M3k[l] + dfl*M3l[l] + force3; 
			fn4[k]  = ff4[k] + df0*M40[l] + df1*M41[l] + df2*M42[l] + df3*M43[l] + df4*M44[l] + df5*M45[l] + df6*M46[l] + dfa*M4a[l] + dfb*M4b[l] + dfc*M4c[l] + dfd*M4d[l] + dfe*M4e[l] + dff*M4f[l] + dfg*M4g[l] + dfh*M4h[l] + dfi*M4i[l] + dfj*M4j[l] + dfk*M4k[l] + dfl*M4l[l] + force4; 
			fn5[k]  = ff5[k] + df0*M50[l] + df1*M51[l] + df2*M52[l] + df3*M53[l] + df4*M54[l] + df5*M55[l] + df6*M56[l] + dfa*M5a[l] + dfb*M5b[l] + dfc*M5c[l] + dfd*M5d[l] + dfe*M5e[l] + dff*M5f[l] + dfg*M5g[l] + dfh*M5h[l] + dfi*M5i[l] + dfj*M5j[l] + dfk*M5k[l] + dfl*M5l[l] + force5; 
			fn6[k]  = ff6[k] + df0*M60[l] + df1*M61[l] + df2*M62[l] + df3*M63[l] + df4*M64[l] + df5*M65[l] + df6*M66[l] + dfa*M6a[l] + dfb*M6b[l] + dfc*M6c[l] + dfd*M6d[l] + dfe*M6e[l] + dff*M6f[l] + dfg*M6g[l] + dfh*M6h[l] + dfi*M6i[l] + dfj*M6j[l] + dfk*M6k[l] + dfl*M6l[l] + force6; 
			fna[k]  = ffa[k] + df0*Ma0[l] + df1*Ma1[l] + df2*Ma2[l] + df3*Ma3[l] + df4*Ma4[l] + df5*Ma5[l] + df6*Ma6[l] + dfa*Maa[l] + dfb*Mab[l] + dfc*Mac[l] + dfd*Mad[l] + dfe*Mae[l] + dff*Maf[l] + dfg*Mag[l] + dfh*Mah[l] + dfi*Mai[l] + dfj*Maj[l] + dfk*Mak[l] + dfl*Mal[l] + forcea; 
			fnb[k]  = ffb[k] + df0*Mb0[l] + df1*Mb1[l] + df2*Mb2[l] + df3*Mb3[l] + df4*Mb4[l] + df5*Mb5[l] + df6*Mb6[l] + dfa*Mba[l] + dfb*Mbb[l] + dfc*Mbc[l] + dfd*Mbd[l] + dfe*Mbe[l] + dff*Mbf[l] + dfg*Mbg[l] + dfh*Mbh[l] + dfi*Mbi[l] + dfj*Mbj[l] + dfk*Mbk[l] + dfl*Mbl[l] + forceb; 
			fnc[k]  = ffc[k] + df0*Mc0[l] + df1*Mc1[l] + df2*Mc2[l] + df3*Mc3[l] + df4*Mc4[l] + df5*Mc5[l] + df6*Mc6[l] + dfa*Mca[l] + dfb*Mcb[l] + dfc*Mcc[l] + dfd*Mcd[l] + dfe*Mce[l] + dff*Mcf[l] + dfg*Mcg[l] + dfh*Mch[l] + dfi*Mci[l] + dfj*Mcj[l] + dfk*Mck[l] + dfl*Mcl[l] + forcec; 
			fnd[k]  = ffd[k] + df0*Md0[l] + df1*Md1[l] + df2*Md2[l] + df3*Md3[l] + df4*Md4[l] + df5*Md5[l] + df6*Md6[l] + dfa*Mda[l] + dfb*Mdb[l] + dfc*Mdc[l] + dfd*Mdd[l] + dfe*Mde[l] + dff*Mdf[l] + dfg*Mdg[l] + dfh*Mdh[l] + dfi*Mdi[l] + dfj*Mdj[l] + dfk*Mdk[l] + dfl*Mdl[l] + forced; 
			fne[k]  = ffe[k] + df0*Me0[l] + df1*Me1[l] + df2*Me2[l] + df3*Me3[l] + df4*Me4[l] + df5*Me5[l] + df6*Me6[l] + dfa*Mea[l] + dfb*Meb[l] + dfc*Mec[l] + dfd*Med[l] + dfe*Mee[l] + dff*Mef[l] + dfg*Meg[l] + dfh*Meh[l] + dfi*Mei[l] + dfj*Mej[l] + dfk*Mek[l] + dfl*Mel[l] + forcee; 
			fnf[k]  = fff[k] + df0*Mf0[l] + df1*Mf1[l] + df2*Mf2[l] + df3*Mf3[l] + df4*Mf4[l] + df5*Mf5[l] + df6*Mf6[l] + dfa*Mfa[l] + dfb*Mfb[l] + dfc*Mfc[l] + dfd*Mfd[l] + dfe*Mfe[l] + dff*Mff[l] + dfg*Mfg[l] + dfh*Mfh[l] + dfi*Mfi[l] + dfj*Mfj[l] + dfk*Mfk[l] + dfl*Mfl[l] + forcef; 
			fng[k]  = ffg[k] + df0*Mg0[l] + df1*Mg1[l] + df2*Mg2[l] + df3*Mg3[l] + df4*Mg4[l] + df5*Mg5[l] + df6*Mg6[l] + dfa*Mga[l] + dfb*Mgb[l] + dfc*Mgc[l] + dfd*Mgd[l] + dfe*Mge[l] + dff*Mgf[l] + dfg*Mgg[l] + dfh*Mgh[l] + dfi*Mgi[l] + dfj*Mgj[l] + dfk*Mgk[l] + dfl*Mgl[l] + forceg; 
			fnh[k]  = ffh[k] + df0*Mh0[l] + df1*Mh1[l] + df2*Mh2[l] + df3*Mh3[l] + df4*Mh4[l] + df5*Mh5[l] + df6*Mh6[l] + dfa*Mha[l] + dfb*Mhb[l] + dfc*Mhc[l] + dfd*Mhd[l] + dfe*Mhe[l] + dff*Mhf[l] + dfg*Mhg[l] + dfh*Mhh[l] + dfi*Mhi[l] + dfj*Mhj[l] + dfk*Mhk[l] + dfl*Mhl[l] + forceh; 
			fni[k]  = ffi[k] + df0*Mi0[l] + df1*Mi1[l] + df2*Mi2[l] + df3*Mi3[l] + df4*Mi4[l] + df5*Mi5[l] + df6*Mi6[l] + dfa*Mia[l] + dfb*Mib[l] + dfc*Mic[l] + dfd*Mid[l] + dfe*Mie[l] + dff*Mif[l] + dfg*Mig[l] + dfh*Mih[l] + dfi*Mii[l] + dfj*Mij[l] + dfk*Mik[l] + dfl*Mil[l] + forcei; 
			fnj[k]  = ffj[k] + df0*Mj0[l] + df1*Mj1[l] + df2*Mj2[l] + df3*Mj3[l] + df4*Mj4[l] + df5*Mj5[l] + df6*Mj6[l] + dfa*Mja[l] + dfb*Mjb[l] + dfc*Mjc[l] + dfd*Mjd[l] + dfe*Mje[l] + dff*Mjf[l] + dfg*Mjg[l] + dfh*Mjh[l] + dfi*Mji[l] + dfj*Mjj[l] + dfk*Mjk[l] + dfl*Mjl[l] + forcej; 
			fnk[k]  = ffk[k] + df0*Mk0[l] + df1*Mk1[l] + df2*Mk2[l] + df3*Mk3[l] + df4*Mk4[l] + df5*Mk5[l] + df6*Mk6[l] + dfa*Mka[l] + dfb*Mkb[l] + dfc*Mkc[l] + dfd*Mkd[l] + dfe*Mke[l] + dff*Mkf[l] + dfg*Mkg[l] + dfh*Mkh[l] + dfi*Mki[l] + dfj*Mkj[l] + dfk*Mkk[l] + dfl*Mkl[l] + forcek; 
			fnl[k]  = ffl[k] + df0*Ml0[l] + df1*Ml1[l] + df2*Ml2[l] + df3*Ml3[l] + df4*Ml4[l] + df5*Ml5[l] + df6*Ml6[l] + dfa*Mla[l] + dfb*Mlb[l] + dfc*Mlc[l] + dfd*Mld[l] + dfe*Mle[l] + dff*Mlf[l] + dfg*Mlg[l] + dfh*Mlh[l] + dfi*Mli[l] + dfj*Mlj[l] + dfk*Mlk[l] + dfl*Mll[l] + forcel; 
    
			gn0[k] = gg0[k] + (ge0 - gg0[k])/tau1;
			gn1[k] = gg1[k] + (ge1 - gg1[k])/tau1;
			gn2[k] = gg2[k] + (ge2 - gg2[k])/tau1;
			gn3[k] = gg3[k] + (ge3 - gg3[k])/tau1;
			gn4[k] = gg4[k] + (ge4 - gg4[k])/tau1;
			gn5[k] = gg5[k] + (ge5 - gg5[k])/tau1;
			gn6[k] = gg6[k] + (ge6 - gg6[k])/tau1;
			gna[k] = gga[k] + (gea - gga[k])/tau1;
			gnb[k] = ggb[k] + (geb - ggb[k])/tau1;
			gnc[k] = ggc[k] + (gec - ggc[k])/tau1;
			gnd[k] = ggd[k] + (ged - ggd[k])/tau1;
			gne[k] = gge[k] + (gee - gge[k])/tau1;
			gnf[k] = ggf[k] + (gef - ggf[k])/tau1;
			gng[k] = ggg[k] + (geg - ggg[k])/tau1;
			gnh[k] = ggh[k] + (geh - ggh[k])/tau1;
			gni[k] = ggi[k] + (gei - ggi[k])/tau1;
			gnj[k] = ggj[k] + (gej - ggj[k])/tau1;
			gnk[k] = ggk[k] + (gek - ggk[k])/tau1;
			gnl[k] = ggl[k] + (gel - ggl[k])/tau1;
		}		  
	}

	//if(t%infoStep==0)
	//	cout << "Process " << rank << ": collided." << endl;
}

