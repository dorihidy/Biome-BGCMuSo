/*
state_init.c
Initialize water, carbon, and nitrogen state variables for pointbgc simulation  

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGCMuSo v7.0.
Original code: Copyright 2000, Peter E. Thornton
Numerical Terradynamic Simulation Group, The University of Montana, USA
Modified code: Copyright 2022, D. Hidy [dori.hidy@gmail.com]
Hungarian Academy of Sciences, Hungary
See the website of Biome-BGCMuSo at http://nimbus.elte.hu/bbgc/ for documentation, model executable and example input files.
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "ini.h"
#include "bgc_struct.h"
#include "bgc_constants.h"
#include "pointbgc_struct.h"
#include "pointbgc_func.h"

int wstate_init(file init, control_struct* ctrl, const siteconst_struct* sitec, const soilprop_struct* sprop, wstate_struct* ws)
{
	int errorCode=0;
	int layer;
	char key[] = "W_STATE";
	char keyword[STRINGSIZE];

	
	/* read water state variable initialization values from *.init */
	if (!errorCode && scan_value(init, keyword, 's'))
	{
		printf("ERROR reading keyword, wstate_init.c\n");
		errorCode=212;
	}
	if (!errorCode && strcmp(keyword,key))
	{
		printf("Expexting keyword --> %s in %s\n",key,init.name);
		errorCode=212;
	}

	if (!errorCode && scan_value(init, &ws->snoww, 'd'))
	{
		printf("ERROR reading snowpack, wstate_init.c\n");
		errorCode=21201;
	}
	if (!errorCode && scan_value(init, &ctrl->wstate_propFC, 'd'))
	{
		printf("ERROR reading initial soilwater (FCprop), wstate_init.c\n");
		errorCode=21202;
	}
	
	/* check that prop_fc is an acceptable proportion  */
	if (!errorCode && (ctrl->wstate_propFC < 0.0))
	{
		printf("ERROR in state_init.c: initial soil water proportion must be >= 0.0 and <= 1.0\n");
		errorCode=2120201;
	}
	if (!errorCode)
	{
		/* calculate initial soilwater in kg/m2 from proportion of
		field capacity volumetric water content, depth, and density of water */
		for (layer = 0; layer < N_SOILLAYERS; layer ++)
		{
			if (ctrl->wstate_propFC > sprop->VWCsat[layer]/sprop->VWCfc[layer])
			{
				printf("ERROR in state_init.c: initial soil water proportion must less than saturation proportion: %lf\n", sprop->VWCsat[layer]/sprop->VWCfc[layer]);
				errorCode=21202;
			}
			else
			{
				ws->soilw[layer] = ctrl->wstate_propFC * sprop->VWCfc[layer] * sitec->soillayer_thickness[layer] * water_density;
				ws->soilw_SUM += ws->soilw[layer];
			}
				
			

		}
	
	}
	
	return (errorCode);
}

int cnstate_init(file init, const epconst_struct* epc, const soilprop_struct* sprop, const siteconst_struct* sitec, 
	             cstate_struct* cs, cinit_struct* cinit, nstate_struct* ns)
{
	int errorCode=0;
	int layer, scanflag, pp;
	int alloc_softstem, alloc_yield, alloc_livestem, alloc_livecroot; 
	char key1[] = "CN_STATE";
	char keyword[STRINGSIZE];
	double NH4_ppm[N_SOILLAYERS];
	double NO3_ppm[N_SOILLAYERS];

	alloc_softstem=alloc_yield=alloc_livestem=alloc_livecroot = 0;

	/* 1. read carbon state variable initial values from *.init */
	if (!errorCode && scan_value(init, keyword, 's'))
	{
		printf("ERROR reading keyword, cstate_init.c\n");
		errorCode=213;
	}
	if (!errorCode && strcmp(keyword,key1))
	{
		printf("Expecting keyword --> %s in %s\n",key1,init.name);
		errorCode=213;
	}
	if (!errorCode && scan_value(init, &cinit->max_leafc, 'd'))
	{
		printf("ERROR reading first-year maximum leaf carbon, cstate_init.c\n");
		errorCode=21301;
	}
	if (!errorCode && scan_value(init, &cinit->max_frootc, 'd'))
	{
		printf("ERROR reading first-year maximum fine root carbon, cstate_init.c\n");
		errorCode=21302;
	}
	if (!errorCode && scan_value(init, &cinit->max_yieldc, 'd'))
	{
		printf("ERROR reading first-year max_yieldc, cstate_init.c\n");
		errorCode=21303;
	}
	if (!errorCode && scan_value(init, &cinit->max_softstemc, 'd'))
	{
		printf("ERROR reading first-year max_sofstemc, cstate_init.c\n");
		errorCode=21304;
	}
	if (!errorCode && scan_value(init, &cinit->max_livestemc, 'd'))
	{
		printf("ERROR reading first-year max_livestemc, cstate_init.c\n");
		errorCode=21305;
	}
	if (!errorCode && scan_value(init, &cinit->max_livecrootc, 'd'))
	{
		printf("ERROR reading first-year max_livecrootc, cstate_init.c\n");
		errorCode=21306;
	}

	/* control */
	for (pp=0; pp<N_PHENPHASES; pp++)
	{
		if (epc->alloc_yield[pp] > 0)     alloc_yield=1;
		if (epc->alloc_softstemc[pp] > 0)  alloc_softstem=1;
		if (epc->alloc_livestemc[pp] > 0)  alloc_livestem=1;
		if (epc->alloc_livecrootc[pp] > 0) alloc_livecroot=1;
	}

	if (alloc_yield == 0)     cinit->max_yieldc = 0;
	if (alloc_softstem == 0)  cinit->max_softstemc = 0;
	if (alloc_livestem == 0)  cinit->max_livestemc = 0;
	if (alloc_livecroot == 0) cinit->max_livecrootc = 0;

	/*--------------------------------------------------*/
	/* 2. read the cwdc initial values in multilayer soil  */

	scanflag=0; 
	for (layer=0; layer<N_SOILLAYERS; layer++)
	{
		if (layer==N_SOILLAYERS-1) scanflag=1;
		if (!errorCode && scan_array(init, &(cs->cwdc[layer]), 'd', scanflag, 1))
		{
			printf("ERROR reading cwdc in layer %i, cstate_init.c\n", layer);
			errorCode=21307;
		}
	}

    /* to avoid dividing by 0: if no deadwood, cwdn is zero. */
	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		if (!errorCode && epc->deadwood_cn > 0.0)
			ns->cwdn[layer] = 0;
		else
			ns->cwdn[layer] = 0;
	}


	/*--------------------------------------------------*/
	/* 3. read the litter carbon pool initial values in multilayer soil  */
	scanflag=0; 
	for (layer=0; layer<N_SOILLAYERS; layer++)
	{
		if (layer==N_SOILLAYERS-1) scanflag=1;
		if (!errorCode && scan_array(init, &(cs->litr1c[layer]), 'd', scanflag, 1))
		{
			printf("ERROR reading litter carbon in labile pool in layer %i, cstate_init.c\n", layer);
			errorCode=21308;
		}
	}
	
	scanflag=0; 
	for (layer=0; layer<N_SOILLAYERS; layer++)
	{
		if (layer==N_SOILLAYERS-1) scanflag=1;
		if (!errorCode && scan_array(init, &(cs->litr2c[layer]), 'd', scanflag, 1))
		{
			printf("ERROR reading litter carbon in unshielded cellulose pool in layer %i, cstate_init.c\n", layer);
			errorCode=21309;
		}
	}
	
	scanflag=0; 
	for (layer=0; layer<N_SOILLAYERS; layer++)
	{
		if (layer==N_SOILLAYERS-1) scanflag=1;
		if (!errorCode && scan_array(init, &(cs->litr3c[layer]), 'd', scanflag, 1))
		{
			printf("ERROR reading litter carbon in shielded cellulose pool in layer %i, cstate_init.c\n", layer);
			errorCode=21310;
		}
	}

	scanflag=0; 
	for (layer=0; layer<N_SOILLAYERS; layer++)
	{
		if (layer==N_SOILLAYERS-1) scanflag=1;
		if (!errorCode && scan_array(init, &(cs->litr4c[layer]), 'd', scanflag, 1))
		{
			printf("ERROR reading litter carbon in lignin pool in layer %i, cstate_init.c\n", layer);
			errorCode=21311;
		}
	}


	/* calculate the litter nitrogen pool initial values for cellulose and and lignin pools, 
	using the leaf litter C:N as the basis for determining N content in all litter components - except of litr1n (reading below)  */
	
	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		ns->litr2n[layer] = cs->litr2c[layer] / epc->leaflitr_cn;
		ns->litr3n[layer] = cs->litr3c[layer] / epc->leaflitr_cn;
		ns->litr4n[layer] = cs->litr4c[layer] / epc->leaflitr_cn;
	}


	/*--------------------------------------------------*/
	/* 4. read the soil carbon pool initial values in multilayer soil  */

	scanflag=0; 
	for (layer=0; layer<N_SOILLAYERS; layer++)
	{
		if (layer==N_SOILLAYERS-1) scanflag=1;
		if (!errorCode && scan_array(init, &(cs->soil1c[layer]), 'd', scanflag, 1))
		{
			printf("ERROR reading labile SOM carbon pool in layer %i, cstate_init.c\n", layer);
			errorCode=21312;
		}
	}
		
	scanflag=0; 
	for (layer=0; layer<N_SOILLAYERS; layer++)
	{
		if (layer==N_SOILLAYERS-1) scanflag=1;
		if (!errorCode && scan_array(init, &(cs->soil2c[layer]), 'd', scanflag, 1))
		{
			printf("ERROR reading fast decomposing SOM carbon pool in layer %i, cstate_init.c\n", layer);
			errorCode=21313;
		}
	}

 	scanflag=0; 
	for (layer=0; layer<N_SOILLAYERS; layer++)
	{
		if (layer==N_SOILLAYERS-1) scanflag=1;
		if (!errorCode && scan_array(init, &(cs->soil3c[layer]), 'd', scanflag, 1))
		{
			printf("ERROR reading slow decomposing SOM carbon pool in layer %i, cstate_init.c\n", layer);
			errorCode=21314;
		}
	}


	scanflag=0; 
	for (layer=0; layer<N_SOILLAYERS; layer++)
	{
		if (layer==N_SOILLAYERS-1) scanflag=1;
		if (!errorCode && scan_array(init, &(cs->soil4c[layer]), 'd', scanflag, 1))
		{
			printf("ERROR reading stable SOM carbon pool in layer %i, cstate_init.c\n", layer);
			errorCode=21315;
		}
	}


	/* multilayer soil */
	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		ns->soil1n[layer] = cs->soil1c[layer]/sprop->soil1_CN;
		ns->soil2n[layer] = cs->soil2c[layer]/sprop->soil2_CN;
		ns->soil3n[layer] = cs->soil3c[layer]/sprop->soil3_CN;
		ns->soil4n[layer] = cs->soil4c[layer]/sprop->soil4_CN;
	}
	
	/* 5. read nitrogen state variable initial values from *.init */
	scanflag=0; 
	for (layer=0; layer<N_SOILLAYERS; layer++)
	{
		if (layer==N_SOILLAYERS-1) scanflag=1;
		if (!errorCode && scan_array(init, &(ns->litr1n[layer]), 'd', scanflag, 1))
		{
			printf("ERROR reading litter nitrogen in labile pool layer %i, cnstate_init.c\n", layer);
			errorCode=21316;
		}
	}
	
	scanflag=0; 
	for (layer=0; layer<N_SOILLAYERS; layer++)
	{
		if (layer==N_SOILLAYERS-1) scanflag=1;
		if (!errorCode && scan_array(init, &(NH4_ppm[layer]), 'd', scanflag, 1))
		{
			printf("ERROR reading soil mineral nitrogen (NH4 pool) in layer %i, cnstate_init.c\n", layer);
			errorCode=21317;
		}
		ns->NH4[layer] = (NH4_ppm[layer] / multi_ppm) * (sprop->BD[layer] * sitec->soillayer_thickness[layer]);
	}

	scanflag=0; 
	for (layer=0; layer<N_SOILLAYERS; layer++)
	{
		if (layer==N_SOILLAYERS-1) scanflag=1;
		if (!errorCode && scan_array(init, &(NO3_ppm[layer]), 'd', scanflag, 1))
		{
			printf("ERROR reading soil mineral nitrogen (NO3 pool) in layer %i, cnstate_init.c\n", layer);
			errorCode=21318;
		}
		ns->NO3[layer] = (NO3_ppm[layer] / multi_ppm) * (sprop->BD[layer] * sitec->soillayer_thickness[layer]);

	}


		
	return (errorCode);
}
