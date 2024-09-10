/* 
groundwater_leaching.c
Calculating soil mineral nitrogen and DOC-DON leaching in groundwater layers 
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGCMuSo v7.0.
Copyright 2022, D. Hidy [dori.hidy@gmail.com]
Hungarian Academy of Sciences, Hungary
See the website of Biome-BGCMuSo at http://nimbus.elte.hu/bbgc/ for documentation, model executable and example input files.
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include "ini.h"
#include "bgc_struct.h"
#include "bgc_func.h"
#include "bgc_constants.h"

int groundwater_leaching(control_struct* ctrl, siteconst_struct* sitec, soilprop_struct* sprop, soilInfo_struct* soilInfo,
	                     cstate_struct* cs, cflux_struct* cf, nstate_struct* ns, nflux_struct* nf, wstate_struct* ws, wflux_struct* wf)
{
	int errorCode, GWlayer, dm;

	double wflux, wstate0, wstate1;
	double state0[N_DISSOLVMATER], state1[N_DISSOLVMATER], conc0[N_DISSOLVMATER], conc1[N_DISSOLVMATER], leachFlux[N_DISSOLVMATER];

	GWlayer = (int)sprop->GWlayer;
	errorCode = 0;


	/*---------------------------------------------------------------------------------*/
	/* 1. Leaching noGW-layer vs zone_NORM/zone_CAPIL*/

	/* 1.1 If there is zone_NORM */
	if (sprop->dz_zoneNORM)
	{
	
		/* 1.1.1 last noGW-layer and zone_NORM*/
		wflux   = wf->soilwFlux[GWlayer - 1];
		
		if (wflux != 0)
		{
			wstate0 = ws->soilw[GWlayer - 1];
			for (dm = 0; dm < N_DISSOLVMATER; dm++) state0[dm] = soilInfo->content_soil[dm][GWlayer - 1];

			wstate1 = sprop->soilw_zoneNORM;
			for (dm = 0; dm < N_DISSOLVMATER; dm++) state1[dm] = soilInfo->content_zoneNORM[dm];


			if (!errorCode && leachCalc(ctrl, soilInfo, wflux, wstate0, wstate1, state0, state1, conc0, conc1, leachFlux))
			{
				printf("\n");
				printf("ERROR in leachsCalc.c for groundwater_leaching.c\n");
				errorCode = 1;
			}

			for (dm = 0; dm < N_DISSOLVMATER; dm++)
			{
				soilInfo->content_soil[dm][GWlayer - 1] = state0[dm];
				soilInfo->conc_soil[dm][GWlayer - 1] = conc0[dm];
				soilInfo->content_zoneNORM[dm] = state1[dm];

				soilInfo->dismatLeach[dm][GWlayer - 1] = leachFlux[dm];

				soilInfo->content_soil[dm][GWlayer] += leachFlux[dm];
			}
		}

		/* 1.1.2 zone_NORM and zone_CAPIL */
		wflux = wf->soilwFlux_NORMvsCAPIL;

		if (wflux != 0)
		{
			wstate0 = sprop->soilw_zoneNORM;
			for (dm = 0; dm < N_DISSOLVMATER; dm++) state0[dm] = soilInfo->content_zoneNORM[dm];

			wstate1 = sprop->soilw_zoneCAPIL;
			for (dm = 0; dm < N_DISSOLVMATER; dm++) state1[dm] = soilInfo->content_zoneCAPIL[dm];


			if (!errorCode && leachCalc(ctrl, soilInfo, wflux, wstate0, wstate1, state0, state1, conc0, conc1, leachFlux))
			{
				printf("\n");
				printf("ERROR in leachsCalc.c for groundwater_leaching.c\n");
				errorCode = 1;
			}

			for (dm = 0; dm < N_DISSOLVMATER; dm++)
			{
				soilInfo->content_zoneNORM[dm] = state0[dm];
				soilInfo->content_zoneCAPIL[dm] = state1[dm];
				soilInfo->dismatLeach_NORM[dm] = leachFlux[dm];
			}
		}

	}
	/* 2.1 If there is only zone_CAPIL */
	else
	{
		wflux = wf->soilwFlux[GWlayer - 1];

		if (wflux != 0)
		{
			wstate0 = ws->soilw[GWlayer - 1];
			for (dm = 0; dm < N_DISSOLVMATER; dm++) state0[dm] = soilInfo->content_soil[dm][GWlayer - 1];

			wstate1 = sprop->soilw_zoneCAPIL;
			for (dm = 0; dm < N_DISSOLVMATER; dm++) state1[dm] = soilInfo->content_zoneCAPIL[dm];


			if (!errorCode && leachCalc(ctrl, soilInfo, wflux, wstate0, wstate1, state0, state1, conc0, conc1, leachFlux))
			{
				printf("\n");
				printf("ERROR in leachCalc.c for groundwater_leaching.c\n");
				errorCode = 1;
			}

			for (dm = 0; dm < N_DISSOLVMATER; dm++)
			{
				soilInfo->content_soil[dm][GWlayer - 1] = state0[dm];
				soilInfo->conc_soil[dm][GWlayer - 1] = conc0[dm];
				soilInfo->content_zoneCAPIL[dm] = state1[dm];

				soilInfo->dismatLeach[dm][GWlayer - 1] = leachFlux[dm];

				soilInfo->content_soil[dm][GWlayer] += leachFlux[dm];
			}
		}

	}

	/*---------------------------------------------------------------------------------*/
	/* 2. Leaching zone_CAPIL vs zone_SAT */

	wflux = wf->GWrecharge - wf->GWdischarge;

	if (wflux != 0)
	{
		wstate0 = sprop->soilw_zoneCAPIL;
		for (dm = 0; dm < N_DISSOLVMATER; dm++) state0[dm] = soilInfo->content_zoneCAPIL[dm];

		wstate1 = sprop->soilw_zoneSAT;
		for (dm = 0; dm < N_DISSOLVMATER; dm++) state1[dm] = soilInfo->content_zoneSAT[dm];


		if (!errorCode && leachCalc(ctrl, soilInfo, wflux, wstate0, wstate1, state0, state1, conc0, conc1, leachFlux))
		{
			printf("\n");
			printf("ERROR in leachsCalc.c for groundwater_leaching.c\n");
			errorCode = 1;
		}

		for (dm = 0; dm < N_DISSOLVMATER; dm++)
		{
			soilInfo->content_zoneCAPIL[dm] = state0[dm];

			if (leachFlux[dm] > 0)
				soilInfo->dismatGWrecharge[dm] = leachFlux[dm];
			else
				soilInfo->dismatGWdischarge[dm] = -1 * leachFlux[dm];

			soilInfo->content_soil[dm][GWlayer] -= leachFlux[dm];
		}
	}

	/* calculation of src/snk variables   */

	for (dm = 0; dm < N_DISSOLVN; dm++)
	{
		ns->GWsnk_N += soilInfo->dismatGWrecharge[dm];
		ns->GWsrc_N += soilInfo->dismatGWdischarge[dm];
	}
	for (dm = N_DISSOLVN; dm < N_DISSOLVMATER; dm++)
	{
		cs->GWsnk_C += soilInfo->dismatGWrecharge[dm];
		cs->GWsrc_C += soilInfo->dismatGWdischarge[dm];
	}

	/* transfer value: content_array ->NH4, NO3, DOC, DON  etc.*/
	if (!errorCode && calc_soilconc(-1, 1, sprop, ws, cs, ns, soilInfo))
	{
		printf("ERROR in calc_soilconc.c for groundwater_leaching.c\n");
		errorCode = 1;
	}

	return (errorCode);
}

