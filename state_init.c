/*
state_init.c
Initialize water, carbon, and nitrogen state variables for pointbgc simulation  

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGCMuSo v7.0.
Original code: Copyright 2000, Peter E. Thornton
Numerical Terradynamic Simulation Group, The University of Montana, USA
Modified code: Copyright 2025, D. Hidy [dori.hidy@gmail.com]
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

int cnstate_init(file init, const control_struct* ctrl, const epconst_struct* epc, const soilprop_struct* sprop, const siteconst_struct* sitec,
	             cstate_struct* cs, cinit_struct* cinit, nstate_struct* ns)
{
	int errorCode=0;
	int layer, scanflag, pp;
	int alloc_softstem, alloc_yield, alloc_livestem, alloc_livecroot; 
	char key1[] = "CN_STATE";
	char keyword[STRINGSIZE];
	double cwdC_ppm[N_SOILLAYERS], cwdN_ppm[N_SOILLAYERS];
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


	scanflag = 0;
	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		if (layer == N_SOILLAYERS - 1) scanflag = 1;
		if (!errorCode && scan_array(init, &(cwdC_ppm[layer]), 'd', scanflag, 1))
		{
			printf("ERROR reading coarse woody debris C content in layer %i, cstate_init.c\n", layer);
			errorCode = 21307;
		}
		cs->cwdc[layer] = (cwdC_ppm[layer] / multi_ppm) * (sprop->BD[layer] * sitec->soillayer_thickness[layer]);
	}



	scanflag = 0;
	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		if (layer == N_SOILLAYERS - 1) scanflag = 1;
		if (!errorCode && scan_array(init, &(cwdN_ppm[layer]), 'd', scanflag, 1))
		{
			printf("ERROR reading coarse woody debris N content in layer %i, cstate_init.c\n", layer);
			errorCode = 21308;
		}
		ns->cwdn[layer] = (cwdN_ppm[layer] / multi_ppm) * (sprop->BD[layer] * sitec->soillayer_thickness[layer]);
	}


	/*--------------------------------------------------*/
	/* 3. read the litter carbon pool initial values in multilayer soil  */
	scanflag = 0;
	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		if (layer == N_SOILLAYERS - 1) scanflag = 1;
		if (!errorCode && scan_array(init, &(cs->litrC_ppm[layer]), 'd', scanflag, 1))
		{
			printf("ERROR reading litter carbon in layer %i, cstate_init.c\n", layer);
			errorCode = 21309;
		}

	}


	scanflag = 0;
	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		if (layer == N_SOILLAYERS - 1) scanflag = 1;
		if (!errorCode && scan_array(init, &(ns->litrN_ppm[layer]), 'd', scanflag, 1))
		{
			printf("ERROR reading litter nitrogen in layer %i, cstate_init.c\n", layer);
			errorCode = 21310;
		}
	}


	scanflag = 0;
	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		if (layer == N_SOILLAYERS - 1) scanflag = 1;
		if (!errorCode && scan_array(init, &(cs->litr4C_ppm[layer]), 'd', scanflag, 1))
		{
			printf("ERROR reading litter carbon in ligning pool in layer %i, cstate_init.c\n", layer);
			errorCode = 21311;
		}
	}


	scanflag = 0;
	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		if (layer == N_SOILLAYERS - 1) scanflag = 1;
		if (!errorCode && scan_array(init, &(ns->litr4N_ppm[layer]), 'd', scanflag, 1))
		{
			printf("ERROR reading litter carbon in ligning pool in layer %i, cstate_init.c\n", layer);
			errorCode = 21312;
		}
	}


	

	/*--------------------------------------------------*/
	/* 4. read the soil carbon pool initial values in multilayer soil  */

	scanflag = 0;
	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		if (layer == N_SOILLAYERS - 1) scanflag = 1;
		if (!errorCode && scan_array(init, &(cs->soilC_ppm[layer]), 'd', scanflag, 1))
		{
			printf("ERROR reading TOC in layer %i, cstate_init.c\n", layer);
			errorCode = 21313;
		}
	}


	scanflag = 0;
	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		if (layer == N_SOILLAYERS - 1) scanflag = 1;
		if (!errorCode && scan_array(init, &(ns->soilN_ppm[layer]), 'd', scanflag, 1))
		{
			printf("ERROR reading TON in layer %i, cstate_init.c\n", layer);
			errorCode = 21314;
		}
	}


	scanflag = 0;
	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		if (layer == N_SOILLAYERS - 1) scanflag = 1;
		if (!errorCode && scan_array(init, &(cs->soil4C_ppm[layer]), 'd', scanflag, 1))
		{
			printf("ERROR reading slow decomposing SOM carbon pool in layer %i, cstate_init.c\n", layer);
			errorCode = 21315;
		}
	}


	scanflag = 0;
	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		if (layer == N_SOILLAYERS - 1) scanflag = 1;
		if (!errorCode && scan_array(init, &(ns->soil4N_ppm[layer]), 'd', scanflag, 1))
		{
			printf("ERROR reading slow decomposing SOM nitrogen in layer %i, cstate_init.c\n", layer);
			errorCode = 21316;
		}
	}


	scanflag = 0;
	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		if (layer == N_SOILLAYERS - 1) scanflag = 1;
		if (!errorCode && scan_array(init, &(NH4_ppm[layer]), 'd', scanflag, 1))
		{
			printf("ERROR reading soil mineral nitrogen (NH4 pool) in layer %i, cnstate_init.c\n", layer);
			errorCode = 21317;
		}
		ns->NH4[layer] = (NH4_ppm[layer] / multi_ppm) * (sprop->BD[layer] * sitec->soillayer_thickness[layer]);
	}

	scanflag = 0;
	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		if (layer == N_SOILLAYERS - 1) scanflag = 1;
		if (!errorCode && scan_array(init, &(NO3_ppm[layer]), 'd', scanflag, 1))
		{
			printf("ERROR reading soil mineral nitrogen (NO3 pool) in layer %i, cnstate_init.c\n", layer);
			errorCode = 21318;
		}
		ns->NO3[layer] = (NO3_ppm[layer] / multi_ppm) * (sprop->BD[layer] * sitec->soillayer_thickness[layer]);

	}

	/* in case of read_restart=0, litter pools are 0, in case of read_restart=1, litter pools are calculated in restart_io.c */

	if (ctrl->spinup == 1)
	{
		for (layer = 0; layer < N_SOILLAYERS; layer++)
		{
			cs->litr1c[layer] = 0;
			cs->litr2c[layer] = 0;
			cs->litr3c[layer] = 0;
			cs->litr4c[layer] = 0;

			ns->litr1n[layer] = 0;
			ns->litr2n[layer] = 0;
			ns->litr3n[layer] = 0;
			ns->litr4n[layer] = 0;

			ns->cwdn[layer] = 0;
			ns->cwdn[layer] = 0;

			cs->soil1c[layer] = 0;
			cs->soil2c[layer] = 0;
			cs->soil3c[layer] = 0;
			cs->soil4c[layer] = 0;

			ns->soil1n[layer] = 0;
			ns->soil2n[layer] = 0;
			ns->soil3n[layer] = 0;
			ns->soil4n[layer] = 0;

			ns->NO3[layer] = 0;
			ns->NH4[layer] = 0;
		}
	}

		
	return (errorCode);
}
