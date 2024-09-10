/* 
sprop_init.c
read sprop file for pointbgc simulation

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGCMuSo v7.0.
Original code: Copyright 2000, Peter E. Thornton
Numerical Terradynamic Simulation Group, The University of Montana, USA
Modified code: Copyright 2022, D. Hidy [dori.hidy@gmail.com]
Hungarian Academy of Sciences, Hungary
See the website of Biome-BGCMuSo at http://nimbus.elte.hu/bbgc/ for documentation, model executable and example input files.
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
20.03.01 Galina Churkina added variable "sum" substituting  t1+t2+t3 in IF statement,
which gave an error.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "ini.h"
#include "bgc_struct.h"
#include "pointbgc_struct.h"
#include "pointbgc_func.h"
#include "bgc_constants.h"
#include "bgc_func.h"


int sprop_init(file init, soilprop_struct* sprop, control_struct* ctrl)
{
	int errorCode=0;
	int layer, scanflag;
	int dofilecloseSOILPROP = 1;
	file sprop_file;
	char key[] = "SOIL_FILE";
	char keyword[STRINGSIZE];
	char header[STRINGSIZE];

	
	/********************************************************************
	**                                                                 **
	** Begin reading initialization file block starting with keyword:  **
	** SOIL_FILE                                                       ** 
	**                                                                 **
	********************************************************************/
	
	
	/* scan for the SOI file keyword, exit if not next */
	if (!errorCode && scan_value(init, keyword, 's'))
	{
		printf("ERROR reading keyword for control data\n");
		errorCode=208;
		dofilecloseSOILPROP = 0;
	}
	if (!errorCode && strcmp(keyword, key))
	{
		printf("Expecting keyword --> %s in file %s\n",key,init.name);
		errorCode=208;
		dofilecloseSOILPROP = 0;
	}
	/* open simple SOI file  */
	if (!errorCode && scan_open(init,&sprop_file,'r',1)) 
	{
		printf("ERROR opening soil file, sprop_init.c\n");
		dofilecloseSOILPROP = 0;
		errorCode=20800;
	}

	/* first scan sprop keyword to ensure proper *.sprop format */
	if (!errorCode && scan_value(sprop_file, header, 's'))
	{
		printf("ERROR reading header, sprop_init.c\n");
		errorCode=208001;
	}

	/***********************************************************/
	/* dividing line from file */ 
	if (!errorCode && scan_value(sprop_file, header, 's'))
	{
		printf("ERROR reading 1. dividing line, sprop_init.c\n");
		errorCode=208002;
	}

	if (!errorCode && scan_value(sprop_file, header, 's'))
	{
		printf("ERROR reading 1. block title, sprop_init.c\n");
		errorCode=208003;
	}
	 
	/****************************/
	/* SOIL GENERIC PARAMETERS */

	/*  soil depth with control*/
	if (!errorCode && scan_value(sprop_file, &sprop->soildepth, 'd'))
	{
		printf("ERROR reading soildepth: sprop_init.c\n");
		errorCode=208004;
	}
	if (!errorCode && sprop->soildepth <= 0)
	{
		printf("ERROR in sprop data in SOI file: negative or zero soildepth, sprop_init.c\n");
		errorCode=2080041;
	}

	if (!errorCode && scan_value(sprop_file, &sprop->soil1_CN, 'd'))
	{
		printf("ERROR reading SOIL1_CN: sprop_init.c\n");
		errorCode=208005;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->soil2_CN, 'd'))
	{
		printf("ERROR reading SOIL2_CN: sprop_init.c\n");
		errorCode=208006;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->soil3_CN, 'd'))
	{
		printf("ERROR reading SOIL3_CN: sprop_init.c\n");
		errorCode=208007;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->soil4_CN, 'd'))
	{
		printf("ERROR reading SOIL4_CN: sprop_init.c\n");
		errorCode=208008;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->critSOCtotal, 'd'))
	{
		printf("ERROR reading critSOCtotal: sprop_init.c\n");
		errorCode=208009;
	}
	
 	if (!errorCode && scan_value(sprop_file, &sprop->NH4_mobilen_prop, 'd'))
	{
		printf("ERROR reading NH4_mobilen_prop: sprop_init.c\n");
		errorCode=208010;
	}
	
	/* aerodynamic paramter (Wallace and Holwill, 1997) */
	if (!errorCode && scan_value(sprop_file, &sprop->aerodyn_resist, 'd'))
	{
		printf("ERROR reading aerodyn_resist: sprop_init.c\n");
		errorCode=208011;
	}
	
	
	/****************************/
	/* dividing line from file */ 

	if (!errorCode && scan_value(sprop_file, header, 's'))
	{
		printf("ERROR reading 2. dividing line, sprop_init.c\n");
		errorCode=208012;
	}
	
	if (!errorCode && scan_value(sprop_file, header, 's'))
	{
		printf("ERROR reading 2. block title, sprop_init.c\n");
		errorCode=208013;
	}
	

	/****************************/
	/*  DECOMPOSITION, NITRIFICATION AND DENITRIFICATION PARAMETERS  */

	if (!errorCode && scan_value(sprop_file, &sprop->VWCratio_DCcrit1, 'd'))
	{
		printf("ERROR reading VWCratio_DCcrit1, sprop_init.c\n");
		errorCode=208014;
	}
	
    if (!errorCode && scan_value(sprop_file, &sprop->VWCratio_DCcrit2, 'd'))
	{
		printf("ERROR reading VWCratio_DCcrit1, sprop_init.c\n");
		errorCode=208015;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->curvature_DC, 'd'))
	{
		printf("ERROR reading curvature_DC, sprop_init.c\n");
		errorCode=208016;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->Tp1_decomp, 'd'))
	{
		printf("ERROR Tp1_decomp, sprop_init.c\n");
		errorCode=208017;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->Tp2_decomp, 'd'))
	{
		printf("ERROR Tp2_decomp, sprop_init.c\n");
		errorCode=208018;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->Tp3_decomp, 'd'))
	{
		printf("ERROR Tp3_decomp, sprop_init.c\n");
		errorCode=208019;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->Tp4_decomp, 'd'))
	{
		printf("ERROR Tp4_decomp, sprop_init.c\n");
		errorCode=208020;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->Tmin_decomp, 'd'))
	{
		printf("ERROR Tp0_decomp, sprop_init.c\n");
		errorCode=208021;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->efolding_depth, 'd'))
	{
		printf("ERROR reading e-folding depth: sprop_init.c\n");
		errorCode=20802;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->netMiner_to_nitrif, 'd'))
	{
		printf("ERROR reading netMiner_to_nitrif: sprop_init.c\n");
		errorCode=208023;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->maxNitrif_rate, 'd'))
	{
		printf("ERROR reading maxNitrif_rate: sprop_init.c\n");
		errorCode=208024;
	}	
	
	if (!errorCode && scan_value(sprop_file, &sprop->N2Ocoeff_nitrif, 'd'))
	{
		printf("ERROR reading N2Ocoeff_nitrif: sprop_init.c\n");
		errorCode=208025;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->pHp1_nitrif, 'd'))
	{
		printf("ERROR pHp1_nitrif, sprop_init.c\n");
		errorCode=208026;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->pHp2_nitrif, 'd'))
	{
		printf("ERROR pHp2_nitrif, sprop_init.c\n");
		errorCode=208027;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->pHp3_nitrif, 'd'))
	{
		printf("ERROR pHp3_nitrif, sprop_init.c\n");
		errorCode=208028;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->pHp4_nitrif, 'd'))
	{
		printf("ERROR pHp4_nitrif, sprop_init.c\n");
		errorCode=208029;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->Tp1_nitrif, 'd'))
	{
		printf("ERROR Tp1_nitrif, sprop_init.c\n");
		errorCode=208030;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->Tp2_nitrif, 'd'))
	{
		printf("ERROR Tp2_nitrif, sprop_init.c\n");
		errorCode=208031;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->Tp3_nitrif, 'd'))
	{
		printf("ERROR Tp3_nitrif, sprop_init.c\n");
		errorCode=208032;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->Tp4_nitrif, 'd'))
	{
		printf("ERROR Tp4_nitrif, sprop_init.c\n");
		errorCode=208033;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->minWFPS_nitrif, 'd'))
	{
		printf("ERROR reading WFPS_min: sprop_init.c\n");
		errorCode=208034;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->opt1WFPS_nitrif, 'd'))
	{
		printf("ERROR reading WFPS_opt1: sprop_init.c\n");
		errorCode=208035;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->opt2WFPS_nitrif, 'd'))
	{
		printf("ERROR reading WFPS_opt2: sprop_init.c\n");
		errorCode=208036;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->scalarWFPSmin_nitrif, 'd'))
	{
		printf("ERROR reading scalarWFPSmin_nitrif: sprop_init.c\n");
		errorCode=208037;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->denitr_coeff, 'd'))
	{
		printf("ERROR reading denitr_coeff: sprop_init.c\n");
		errorCode=208038;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->N2Oratio_denitr, 'd'))
	{
		printf("ERROR reading N2Oratio_denitr: sprop_init.c\n");
		errorCode=208039;
	}
	if (!errorCode && scan_value(sprop_file, &sprop->critWFPS_denitr, 'd'))
	{
		printf("ERROR reading critWFPS_denitr: sprop_init.c\n");
		errorCode=208040;
	}
	
	
	/****************************/
	/* dividing line from file */ 

	if (!errorCode && scan_value(sprop_file, header, 's'))
	{
		printf("ERROR reading 3. dividing line, sprop_init.c\n");
		errorCode=208041;
	}
	
	if (!errorCode && scan_value(sprop_file, header, 's'))
	{
		printf("ERROR reading 3. block title, sprop_init.c\n");
		errorCode=208042;
	}
	
	
	/****************************/
	/* RATE SCALARS */

	/* respiration fractions for fluxes between compartments  */
	if (!errorCode && scan_value(sprop_file, &sprop->rfl1s1, 'd'))
	{
		printf("ERROR reading rfl1s1: sprop_init.c\n");
		errorCode=208043;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->rfl2s2, 'd'))
	{
		printf("ERROR reading rfl2s2: sprop_init.c\n");
		errorCode=208044;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->rfl4s3, 'd'))
	{
		printf("ERROR reading rfl4s3: sprop_init.c\n");
		errorCode=208045;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->rfs1s2, 'd'))
	{
		printf("ERROR reading rfs1s2: sprop_init.c\n");
		errorCode=208046;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->rfs2s3, 'd'))
	{
		printf("ERROR reading rfs2s3: sprop_init.c\n");
		errorCode=208047;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->rfs3s4, 'd'))
	{
		printf("ERROR reading rfs3s4: sprop_init.c\n");
		errorCode=208048;
	}
	
	/* 	base values of rate constants are (1/day)   */
	if (!errorCode && scan_value(sprop_file, &sprop->kl1_base, 'd'))
	{
		printf("ERROR reading kl1_base: sprop_init.c\n");
		errorCode=208049;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->kl2_base, 'd'))
	{
		printf("ERROR reading kl2_base: sprop_init.c\n");
		errorCode=208050;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->kl4_base, 'd'))
	{
		printf("ERROR reading kl4_base: sprop_init.c\n");
		errorCode=208051;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->ks1_base, 'd'))
	{
		printf("ERROR reading ks1_base: sprop_init.c\n");
		errorCode=208052;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->ks2_base, 'd'))
	{
		printf("ERROR reading ks2_base: sprop_init.c\n");
		errorCode=208053;
	}
	if (!errorCode && scan_value(sprop_file, &sprop->ks3_base, 'd'))
	{
		printf("ERROR reading ks3_base: sprop_init.c\n");
		errorCode=208054;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->ks4_base, 'd'))
	{
		printf("ERROR reading ks4_base: sprop_init.c\n");
		errorCode=208055;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->kfrag_base, 'd'))
	{
		printf("ERROR reading kfrag_base: sprop_init.c\n");
		errorCode=208056;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->L1release_ratio, 'd'))
	{
		printf("ERROR reading L1release_ratio: sprop_init.c\n");
		errorCode=208057;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->L2release_ratio, 'd'))
	{
		printf("ERROR reading L2release_ratio: sprop_init.c\n");
		errorCode=208058;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->L4release_ratio, 'd'))
	{
		printf("ERROR reading L4release_ratio: sprop_init.c\n");
		errorCode=208059;
	}
	

	/****************************/
	/* dividing line from file */ 

	if (!errorCode && scan_value(sprop_file, header, 's'))
	{
		printf("ERROR reading 4. dividing line, sprop_init.c\n");
		errorCode=208060;
	}
	
	if (!errorCode && scan_value(sprop_file, header, 's'))
	{
		printf("ERROR reading 4. block title, sprop_init.c\n");
		errorCode=208061;
	}
	

	/****************************/
	/* SOIL MOISTURE PARAMETERS */

	/*  ratio of bare soil evaporation and pot.evaporation */
	if (!errorCode && scan_value(sprop_file, &sprop->soilEVPlim, 'd'))
	{
		printf("ERROR reading soilEVPlim: sprop_init.c\n");
		errorCode=208062;
	}
	
	/*  maximum height of pond water */
	if (!errorCode && scan_value(sprop_file, &sprop->pondmax, 'd'))
	{
		printf("ERROR reading maximum height of pond water: sprop_init.c\n");
		errorCode=208063;
	}
	
	/*  curvature of soilstress function */
	if (!errorCode && scan_value(sprop_file, &sprop->curvature_WS, 'd'))
	{
		printf("ERROR reading curvature_WS: sprop_init.c\n");
		errorCode=208064;
	}
	
	/* runoff parameter (Campbell and Diaz)  */
	if (!errorCode && scan_value(sprop_file, &sprop->RCN_mes, 'd'))
	{
		printf("ERROR reading measured runoff curve number: sprop_init.c\n");
		errorCode=208065;
	}
	
	/* fraction of dissolving coefficients  */
	if (!errorCode && scan_value(sprop_file, &sprop->SOIL1_dissolv_prop, 'd'))
	{
		printf("ERROR reading SOIL1_dissolv_prop: sprop_init.c\n");
		errorCode=208066;
	}	
	
	if (!errorCode && scan_value(sprop_file, &sprop->SOIL2_dissolv_prop, 'd'))
	{
		printf("ERROR reading SOIL2_dissolv_prop: sprop_init.c\n");
		errorCode=208067;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->SOIL3_dissolv_prop, 'd'))
	{
		printf("ERROR reading SOIL3_dissolv_prop: sprop_init.c\n");
		errorCode=208068;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->SOIL4_dissolv_prop, 'd'))
	{
		printf("ERROR reading SOIL4_dissolv_prop: sprop_init.c\n");
		errorCode=208069;
	}
	
	/* soilCover parameters  */
	if (!errorCode && scan_value(sprop_file, &sprop->pLAY_soilCover, 'd'))
	{
		printf("ERROR pLAY_soilCover, sprop_init.c\n");
		errorCode=208070;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->pCRIT_soilCover, 'd'))
	{
		printf("ERROR pCRIT_soilCover, sprop_init.c\n");
		errorCode=208071;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->p1_soilCover, 'd'))
	{
		printf("ERROR p1_soilCover, sprop_init.c\n");
		errorCode=208072;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->p2_soilCover, 'd'))
	{
		printf("ERROR p2_soilCover, sprop_init.c\n");
		errorCode=208073;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->p3_soilCover, 'd'))
	{
		printf("ERROR p3_soilCover, sprop_init.c\n");
		errorCode=208074;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->pRED_soilCover, 'd'))
	{
		printf("ERROR pRED_soilCover, sprop_init.c\n");
		errorCode=208075;
	}
	
	/* tipping diffusion parameters  */
	if (!errorCode && scan_value(sprop_file, &sprop->p1diffus_tipping, 'd'))
	{
		printf("ERROR p1diffus_tipping, sprop_init.c\n");
		errorCode=208076;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->p2diffus_tipping, 'd'))
	{
		printf("ERROR p2diffus_tipping, sprop_init.c\n");
		errorCode=208077;
	}
	
	if (!errorCode && scan_value(sprop_file, &sprop->p3diffus_tipping, 'd'))
	{
		printf("ERROR p3diffus_tipping, sprop_init.c\n");
		errorCode=208078;
	}
	
	
	/****************************/
	/* dividing line from file */ 
	if (!errorCode && scan_value(sprop_file, header, 's'))
	{
		printf("ERROR reading 5. dividing line, sprop_init.c\n");
		errorCode=208079;
	}
	
	if (!errorCode && scan_value(sprop_file, header, 's'))
	{
		printf("ERROR reading 5. block title, sprop_init.c\n");
		errorCode=208080;
	}
	

	/****************************/
	/* CH4 PARAMETERS */

    if (!errorCode && scan_value(sprop_file, &sprop->pBD1_CH4, 'd'))
	{
		printf("ERROR reading pBD1_CH4: sprop_init.c\n");
		errorCode=208081;
	}
	

	if (!errorCode && scan_value(sprop_file, &sprop->pBD2_CH4, 'd'))
	{
		printf("ERROR reading pBD2_CH4: sprop_init.c\n");
		errorCode=208082;
	}
	

	if (!errorCode && scan_value(sprop_file, &sprop->pVWC1_CH4, 'd'))
	{
		printf("ERROR reading pVWC1_CH4: sprop_init.c\n");
		errorCode=208083;
	}
	

	if (!errorCode && scan_value(sprop_file, &sprop->pVWC2_CH4, 'd'))
	{
		printf("ERROR reading pVWC2_CH4: sprop_init.c\n");
		errorCode=208084;
	}
	

	if (!errorCode && scan_value(sprop_file, &sprop->pVWC3_CH4, 'd'))
	{
		printf("ERROR reading pVWC3_CH4: sprop_init.c\n");
		errorCode=208085;
	}
	

	if (!errorCode && scan_value(sprop_file, &sprop->C_pVWC4, 'd'))
	{
		printf("ERROR reading C_pVWC4: sprop_init.c\n");
		errorCode=208086;
	}
	

	if (!errorCode && scan_value(sprop_file, &sprop->pTS_CH4, 'd'))
	{
		printf("ERROR reading pTS_CH4: sprop_init.c\n");
		errorCode=208087;
	}
	
	
	/****************************/
	/* dividing line from file */ 
	if (!errorCode && scan_value(sprop_file, header, 's'))
	{
		printf("ERROR reading 6. dividing line, sprop_init.c\n");
		errorCode=208088;
	}
	
	if (!errorCode && scan_value(sprop_file, header, 's'))
	{
		printf("ERROR reading 6. block title, sprop_init.c\n");
		errorCode=208089;
	}
	
	
	/****************************/
	/* SOIL COMPOSITION AND CHARACTERISTIC VALUES (-9999: no measured data) */
	
	/* SAND array - mulilayer soil  */
	scanflag=0; 
	for (layer=0; layer<N_SOILLAYERS; layer++)
	{
		if (layer==N_SOILLAYERS-1) scanflag=1;
		if (!errorCode && scan_array(sprop_file, &(sprop->sand[layer]), 'd', scanflag, 1))
		{
			printf("ERROR reading percent sand in layer %i, sprop_init.c\n", layer);
			errorCode=208090;
		}
	}
	

	/* SILT array - mulilayer soil   */
	scanflag=0; 
	for (layer=0; layer<N_SOILLAYERS; layer++)
	{
		if (layer==N_SOILLAYERS-1) scanflag=1;
		if (!errorCode && scan_array(sprop_file, &(sprop->silt[layer]), 'd', scanflag, 1))
		{
			printf("ERROR reading percent silt in layer %i, sprop_init.c\n", layer);
			errorCode=208091;
		}
	}
	

	/* pH array - mulilayer soil */
	scanflag=0; 
	for (layer=0; layer<N_SOILLAYERS; layer++)
	{
		if (layer==N_SOILLAYERS-1) scanflag=1;
		if (!errorCode && scan_array(sprop_file, &(sprop->pH[layer]), 'd', scanflag, 1))
		{
			printf("ERROR reading soil pH in layer %i, sprop_init.c\n", layer);
			errorCode=208092;
		}
	}
	

	/* soilB parameter - mulilayer soil */
	scanflag=0; 
	for (layer=0; layer<N_SOILLAYERS; layer++)
	{
		if (layer==N_SOILLAYERS-1) scanflag=1;
		if (!errorCode && scan_array(sprop_file, &(sprop->soilB[layer]), 'd', scanflag, 1))
		{
			printf("ERROR reading soilB in layer %i, sprop_init.c\n", layer);
			errorCode=208093;
		}
	}
	

	/* measured bulk density    */
	scanflag=0; 
	for (layer=0; layer<N_SOILLAYERS; layer++)
	{
		if (layer==N_SOILLAYERS-1) scanflag=1;
		if (!errorCode && scan_array(sprop_file, &(sprop->BD_mes[layer]), 'd', scanflag, 1))
		{
			printf("ERROR reading BD_mes in layer %i, sprop_init.c\n", layer);
			errorCode=208094;
		}
	}
	

	/* measured critical VWC values - saturation    */
	scanflag=0; 
	for (layer=0; layer<N_SOILLAYERS; layer++)
	{
		if (layer==N_SOILLAYERS-1) scanflag=1;
		if (!errorCode && scan_array(sprop_file, &(sprop->VWCsat_mes[layer]), 'd', scanflag, 1))
		{
			printf("ERROR reading VWCsat_mes in layer %i, sprop_init.c\n", layer);
			errorCode=208095;
		}
	}
	
	
	/* measured critical VWC values - field capacity     */
	scanflag=0; 
	for (layer=0; layer<N_SOILLAYERS; layer++)
	{
		if (layer==N_SOILLAYERS-1) scanflag=1;
		if (!errorCode && scan_array(sprop_file, &(sprop->VWCfc_mes[layer]), 'd', scanflag, 1))
		{
			printf("ERROR reading VWCfc_mes in layer %i, sprop_init.c\n", layer);
			errorCode=208096;
		}
	}
	
	
	/* measured critical VWC values - wilting point    */
	scanflag=0; 
	for (layer=0; layer<N_SOILLAYERS; layer++)
	{
		if (layer==N_SOILLAYERS-1) scanflag=1;
		if (!errorCode && scan_array(sprop_file, &(sprop->VWCwp_mes[layer]), 'd', scanflag, 1))
		{
			printf("ERROR reading VWCwp_mes in layer %i, sprop_init.c\n", layer);
			errorCode=208097;
		}
	}
	

	/* measured critical VWC values - hygr. water    */
	scanflag=0; 
	for (layer=0; layer<N_SOILLAYERS; layer++)
	{
		if (layer==N_SOILLAYERS-1) scanflag=1;
		if (!errorCode && scan_array(sprop_file, &(sprop->VWChw_mes[layer]), 'd', scanflag, 1))
		{
			printf("ERROR reading VWChw_mes in layer %i, sprop_init.c\n", layer);
			errorCode=208098;
		}
	}
	
	/* CONTROL - measured VWC values (all or none should to be set) */
	if (sprop->VWCsat_mes[0] == DATA_GAP || sprop->VWCfc_mes[0] == DATA_GAP || sprop->VWCwp_mes[0] == DATA_GAP  || sprop->VWChw_mes[0] == DATA_GAP)
	{	
		if (sprop->VWCsat_mes[0] == DATA_GAP && sprop->VWCfc_mes[0] == DATA_GAP && sprop->VWCwp_mes[0] == DATA_GAP  && sprop->VWChw_mes[0] == DATA_GAP)
		{	
			for (layer=0; layer < N_SOILLAYERS; layer++)
			{
				sprop->VWCsat_mes[layer] = DATA_GAP;
				sprop->VWCfc_mes[layer]  = DATA_GAP;
				sprop->VWCwp_mes[layer]  = DATA_GAP;
				sprop->VWChw_mes[layer]  = DATA_GAP;
			}
		}
		else
		{
			printf("ERROR in measured soil data in SOI file\n");
			printf("All or none critical VWC data should to be set by the user\n");
			errorCode=2080981;
		}
	}
	

	/* measured drainage coeff */
	scanflag=0; 
	for (layer=0; layer<N_SOILLAYERS; layer++)
	{
		if (layer==N_SOILLAYERS-1) scanflag=1;
		if (!errorCode && scan_array(sprop_file, &(sprop->drainCoeff_mes[layer]), 'd', scanflag, 1))
		{
			printf("ERROR reading drainCoeff_mes in layer %i, sprop_init.c\n", layer);
			errorCode=208099;
		}
	}
	

	/* measured hydraulic conductivity  */
	scanflag=0; 
	for (layer=0; layer<N_SOILLAYERS; layer++)
	{
		if (layer==N_SOILLAYERS-1) scanflag=1;
		if (!errorCode && scan_array(sprop_file, &(sprop->hydrCONTUCTsatMES_cmPERday[layer]), 'd', scanflag, 1))
		{
			printf("ERROR reading hydrCONTUCTsatMES_cmPERday in layer %i, sprop_init.c\n", layer);
			errorCode=208100;
		}
	}
	

   /* capillary fringe */
	scanflag=0; 
	for (layer=0; layer<N_SOILLAYERS; layer++)
	{
		if (layer==N_SOILLAYERS-1) scanflag=1;
		if (!errorCode && scan_array(sprop_file, &(sprop->CapillFringeMES_cm[layer]), 'd', scanflag, 1))
		{
			printf("ERROR reading CapillFringeMES_cm in layer %i, sprop_init.c\n", layer);
			errorCode=208101;
		}
	}
	
	


	/* MULTILAYER SOIL CALCULATION: calculation of soil properties layer by layer (soilb, BD, PSI and VWC critical values),hydr.conduct and hydr.diffus  */
	if (!errorCode && multilayer_soilcalc(ctrl, sprop))
	{
		printf("\n");
		printf("ERROR in multilayer_soilcalc.c for sprop_init.c\n");
		errorCode=208; 
	} 

	/* groundwater data initalization */
	sprop->GWD = DATA_GAP;
	sprop->FLD = DATA_GAP;
	sprop->GWlayer=DATA_GAP;


	/* -------------------------------------------*/
	if (dofilecloseSOILPROP) fclose(sprop_file.ptr);

		
	return (errorCode);

}

int soilb_estimation(double sand, double silt, double* soilB, double* VWCsat,double* VWCfc, double* VWCwp,  
	                 double* BD, double* RCN, int* soiltype)
{

	int st=-1;
	int errorCode=0;
	
	double soilb_array[12]			= {3.45,  4.11,  5.26,  6.12,  5.39,  4.02,		7.63,  7.71,	8.56,	9.22,	10.45,	12.46};
	double VWCsat_array[12]			= {0.4,   0.42,  0.44,  0.46,  0.48,  0.49,		0.5,   0.505,	0.51,	0.515,	0.52,	0.525};
	double VWCfc_array[12]			= {0.155, 0.190, 0.250, 0.310, 0.360, 0.380,	0.390, 0.405,	0.420,	0.435,	0.445,	0.460};
	double VWCwp_array[12]			= {0.030, 0.050, 0.090, 0.130, 0.170, 0.190,	0.205, 0.220,	0.240,	0.260,	0.275,	0.290};
	double BD_array[12]				= {1.6,   1.58,  1.56,  1.54,  1.52,  1.5,		1.48,  1.46,	1.44,	1.42,	1.4,	1.38};
	double RCN_array[12]			= {50,    52,    54,    56,    58,    60,		62,    64,		66,		68,		70,		72};
	                                 //sand l.sand sa.loam  loam  si.loam silt  sa.c.loam  c.loam si.c.loam sa.clay si.clay clay  
	double clay = 100-sand-silt;



 	if (silt+1.5*clay < 15)													//sand 
		st=0;	

	if ((silt+1.5*clay >= 15) && (silt+2*clay < 30)) 						//loamy_sand
		st=1;	

	if ((clay >= 7 && clay < 20) && (sand > 52) && (silt+2*clay >= 30))		//sandy_loam
		st=2;  

	if (clay < 7 && silt < 50 && silt+2*clay >= 30)   						//sandy_loam
		st=2;

	if (clay >= 7 && clay < 27 && silt >= 28 && silt < 50 && sand <= 52)	//loam
		st=3;

	if (silt >= 50 && clay >= 12 && clay < 27)								//silt loam
		st=4;

	if (silt >= 50 && silt < 80 && clay < 12)								//silt loam
		st=4;

	if (silt >= 80 && clay < 12)											//silt
		st=5;

	if (clay >= 20 && clay < 35 && silt < 28 && sand > 45)					//sandy_clay_loam
		st=6;

	if (clay >= 27 && clay < 40 && sand > 20 && sand <= 45) 				//clay_loam	 
		st=7;

	if (clay >= 27 && clay < 40 && sand <= 20)  							//silty clay_loam				 
		st=8;

	if (clay >= 35 && sand > 45)  											//sandy_clay					 
		st=9;

	if (clay >= 40 && silt >= 40) 											//silty_clay						 
 		st=10;

	if (clay >= 40 && sand <= 45 && silt < 40) 								//clay					 
		st=11;

	*soiltype = st;

	if (st >= 0)
	{
		*soilB			= soilb_array[st];
		*VWCsat			= VWCsat_array[st];
		*VWCfc			= VWCfc_array[st];
		*VWCwp			= VWCwp_array[st];
		*BD				= BD_array[st];
		*RCN			= RCN_array[st];
	}
	else
	{	
		printf("ERROR in soiltype founding\n");
		errorCode=1;
	}

  return (errorCode);
  
}

int multilayer_soilcalc(control_struct* ctrl,  soilprop_struct* sprop)
{
	int layer;

	double sand, silt, clay, PSIfc, PSIwp, PSIsat;
	double soilB, BD, RCN, CapillFringe, VWCsat, VWCfc, VWCwp, VWChw, hydrCONDUCTsat, hydrDIFFUSsat, hydrCONDUCTfc,hydrDIFFUSfc;
	int errorCode = 0;
	int soiltype = 0;

	soilB = BD = RCN = CapillFringe = PSIsat = VWCsat =VWCwp = VWChw = hydrCONDUCTsat = hydrDIFFUSsat = hydrCONDUCTfc =hydrDIFFUSfc = 0;

	
	/* -------------------------------------------------------------------------------------------------------------------------------*/

	

	/* 1. estimated soil water potential at hygroscopic water in MPa (1MPa = 10000cm)  (fc: pF = 2.5; wp: pF = 4.2) */
	sprop->PSIhw  = pow(10,pF_hygroscopw) / (-10000);

	
	/* 2. Calculate the soil pressure-volume coefficients from texture data (soil water content, soil water potential and Clapp-Hornberger parameter) 
		    - different estimation methods 4: modell estimtaion, talajharomszog, measured data )*/

	for (layer=0; layer < N_SOILLAYERS; layer++)
	{

		sand		= sprop->sand[layer];
		silt		= sprop->silt[layer];
		clay		= 100-sand-silt;

		sprop->clay[layer] = clay;

	
		/*  2.1 CONTROL to avoid negative data  */
 		if (sand < 0 || silt < 0 || clay < 0)
		{
			if (!errorCode) printf("ERROR in site data in SOI file: negative sand/silt/clay data, sprop_init.c\n");
			if (!errorCode) printf("Check values in initialization file.\n");
			errorCode=1;
		}

		if (sand == 0 && silt == 0 && sprop->pH[layer] == 0)
		{
			if (!errorCode) printf("ERROR in site data in SOI file: user-defined sand and silt and pH data is zero in layer %i\n", layer);
			if (!errorCode) printf("Check values in initialization file.\n");
			errorCode=1;
		}
	

	
		/* 2.2 saturation value of soil water potential */
		 PSIsat = -(exp((1.54 - 0.0095*sand + 0.0063*silt)*log(10.0))*9.8e-5);
		
		 /* capillary fringe */
		if (sprop->CapillFringeMES_cm[layer] == DATA_GAP) 
			sprop->CapillFringe[layer] = -100*PSIsat;
		else
			sprop->CapillFringe[layer] = sprop->CapillFringeMES_cm[layer] / 100; // dimension from cm to m

		if (soilb_estimation(sand, silt, &soilB, &VWCsat, &VWCfc, &VWCwp,&BD, &RCN, &soiltype))
		{
			if (!errorCode) 
			{
				printf("\n");
				printf("ERROR in soilb_estimation.c for sprop_init\n");
			}
			errorCode=1;
		}	
		ctrl->soiltype_array[layer] = soiltype;
	
		if (sprop->BD_mes[layer] != (double) DATA_GAP) BD      = sprop->BD_mes[layer];

		/* using the model soil properties */
		if (sprop->VWCsat_mes[layer] == (double) DATA_GAP)
		{
	        /* VWCfc, VWCwp VWChw are the functions of VWCsat and soilB (new HW function from Fodor Nándor */
			VWChw= (0.01 * (100 - sand) + 0.06 * clay)/100.;
			//VWChw = VWCsat * (log(soilB) / log(sprop->PSIhw/PSIsat));
		    PSIfc = exp(VWCsat/VWCfc*log(soilB))*PSIsat;
			PSIwp = exp(VWCsat/VWCwp*log(soilB))*PSIsat;


			/* control for soil type with high clay content  */
			if (VWCsat - VWCfc < 0.001)  VWCfc = VWCsat - 0.001;
			if (VWCfc  - VWCwp < 0.001)  VWCwp = VWCfc  - 0.001;
			if (VWCwp  - VWChw < 0.001)  VWChw = VWCwp  - 0.001;
	
		}
		/* measured soil water content at wilting point, field capacity and saturation in m3/m3 */
		else 	
		{

			VWCsat = sprop->VWCsat_mes[layer];
			VWCfc  = sprop->VWCfc_mes[layer];
			VWCwp  = sprop->VWCwp_mes[layer];
			VWChw  = sprop->VWChw_mes[layer];

			PSIfc = exp(VWCsat/VWCfc*log(soilB))*PSIsat;
			PSIwp = exp(VWCsat/VWCwp*log(soilB))*PSIsat;
			
		}
		
	
		/* 2.4 CONTROL - measured VWC values: SAT>FC>WP>HW */
		if ((VWCsat - VWCfc) < 0.0001 || (VWCfc - VWCwp) < 0.0001  || (VWCwp - VWChw) < 0.0001 || VWChw < 0.0001 || VWCsat > 1.0) 
		{
			if (!errorCode) printf("ERROR in measured VWC data in SOI file\n");
			if (!errorCode) printf("rules: VWCsat > VWCfc; VWCfc > VWCwp; VWCwp > VWChw; VWCsat <1.0, VWChw>0.01\n");
			errorCode=1;
		}


			

		if (sprop->soilB[layer] == DATA_GAP) sprop->soilB[layer] = soilB;
		sprop->BD[layer]                = BD;
		sprop->PSIsat[layer]			= PSIsat;
		sprop->PSIfc[layer]				= PSIfc;
		sprop->PSIwp[layer]				= PSIwp;
		sprop->VWCsat[layer]			= VWCsat;
		sprop->VWCfc[layer]				= VWCfc;
		sprop->VWCfc_base[layer]		= VWCfc;
		sprop->VWCwp[layer]				= VWCwp;
		sprop->VWChw[layer]				= VWChw;
		sprop->hydrCONDUCTsat[layer]	= hydrCONDUCTsat;
		sprop->hydrDIFFUSsat[layer]		= hydrDIFFUSsat;
		sprop->hydrCONDUCTfc[layer]		= hydrCONDUCTfc;
		sprop->hydrDIFFUSfc[layer]		= hydrDIFFUSfc;


		if (layer == 0) 
		{
			if (sprop->RCN_mes == DATA_GAP)
				sprop->RCN  = RCN;
			else
				sprop->RCN  = sprop->RCN_mes;
		}

		/* soil water conductitvity constans: ratio of the drained of water of a given day [1/day] */

		
	}

	
	return (errorCode);

}