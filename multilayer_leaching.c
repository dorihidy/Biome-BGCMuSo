/* 
multilayer_leaching.c
Calculating soil mineral nitrogen and DOC-DON leaching in multilayer soil 
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

int multilayer_leaching(control_struct* ctrl, soilprop_struct* sprop, soilInfo_struct* soilInfo, 
	                    cstate_struct* cs, cflux_struct* cf, nstate_struct* ns, nflux_struct* nf, wstate_struct* ws, wflux_struct* wf)
{
	int errorCode, layer,  N_noGWlayers, dm;

	double wflux, wstate0, wstate1;
	double state0[N_DISSOLVMATER], state1[N_DISSOLVMATER], conc0[N_DISSOLVMATER], conc1[N_DISSOLVMATER], leachFlux[N_DISSOLVMATER];

	errorCode = 0; 

	/*---------------------------------------------------------------------------------*/
	/* multilayer_leaching.c is used only for layers without GW */
	if (sprop->GWD == DATA_GAP)
		N_noGWlayers = N_SOILLAYERS;
	else
		N_noGWlayers = (int)sprop->GWlayer;



	/*---------------------------------------------------------------------------------*/
	/* leaching fluxes for the 10 dissolving material types */
	for (layer = 0; layer < N_noGWlayers - 1; layer++)
	{

		wflux = wf->soilwFlux[layer];
		if (wflux != 0)
		{
			wstate0 = ws->soilw[layer];
			for (dm = 0; dm < N_DISSOLVMATER; dm++) state0[dm] = soilInfo->content_soil[dm][layer];
			if (layer + 1 < N_SOILLAYERS)
			{
				wstate1 = ws->soilw[layer + 1];
				for (dm = 0; dm < N_DISSOLVMATER; dm++) state1[dm] = soilInfo->content_soil[dm][layer + 1];
			}
			else
			{
				wstate1 = ws->soilw[layer];
				for (dm = 0; dm < N_DISSOLVMATER; dm++) state1[dm] = soilInfo->content_soil[dm][layer];
			}

			if (!errorCode && leachCalc(ctrl, soilInfo, wflux, wstate0, wstate1, state0, state1, conc0, conc1, leachFlux))
			{
				printf("\n");
				printf("ERROR in leachCalc.c for muliltayer_leaching.c\n");
				errorCode = 1;
			}

			for (dm = 0; dm < N_DISSOLVMATER; dm++)
			{
				soilInfo->content_soil[dm][layer] = state0[dm];
				soilInfo->conc_soil[dm][layer] = conc0[dm];
				if (layer + 1 < N_SOILLAYERS)
				{
					soilInfo->content_soil[dm][layer + 1] = state1[dm];
					soilInfo->conc_soil[dm][layer + 1] = conc1[dm];
				}
				soilInfo->dismatLeach[dm][layer] = leachFlux[dm];
			}


		}
	}

	/*---------------------------------------------------------------------------------*/
	/* transfer value: content_array ->NH4, NO3, DOC, DON  etc.*/
	if (!errorCode && calc_soilconc(-1, 1, sprop, ws, cs, ns, soilInfo))
	{
		printf("ERROR in calc_soilconc.c for multilayer_leaching.c\n");
		errorCode = 1;
	}


	/*---------------------------------------------------------------------------------*/
	/* deepleach calculation from the bottom layer */

	for (dm = 0; dm < N_DISSOLVN; dm++)              ns->Ndeepleach_snk += soilInfo->dismatLeach[dm][N_SOILLAYERS - 1];
	for (dm = N_DISSOLVN; dm < N_DISSOLVMATER; dm++) cs->Cdeepleach_snk += soilInfo->dismatLeach[dm][N_SOILLAYERS - 1];



	
	
	return (errorCode);
}

