/* 
multilayer_leaching.c
Calculating soil mineral nitrogen and DOC-DON leaching in multilayer soil 
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGCMuSo v7.0.
Copyright 2025, D. Hidy [dori.hidy@gmail.com]
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

int multilayer_leaching(const siteconst_struct* sitec, soilprop_struct* sprop, soilInfo_struct* soilInfo, cstate_struct* cs,  nstate_struct* ns, wstate_struct* ws, wflux_struct* wf)
{
	int errorCode, layer, GWlayer, CFlayer, dm;

	double wflux, diffNORM, diffCAPIL, diff;
	double change, dismatLeachNORM, dismatLeach_NORMvsCAPIL, dismatLeach_NORMfromAbove, dismatLeach_CAPILfromAbove, dischargeNORM, dischargeCAPIL, rechargeNORM, rechargeCAPIL, dismatLeachCAPIL, percolDiffus_NORM, dismatLeach_fromAbove;
	double soilwAVAIL[N_DISSOLVMATER], conc0[N_DISSOLVMATER], conc1[N_DISSOLVMATER];

	errorCode = 0; 

	change= dismatLeachNORM = wflux = dismatLeach_NORMvsCAPIL = dismatLeach_NORMfromAbove = dismatLeach_CAPILfromAbove = dischargeNORM = dischargeCAPIL = rechargeNORM = rechargeCAPIL = dismatLeachCAPIL = percolDiffus_NORM = dismatLeach_fromAbove = 0;
	GWlayer = (int)sprop->GWlayer;
	CFlayer = (int)sprop->CFlayer;

	/*---------------------------------------------------------------------------------*/
	/* for concentration calculation original soilw data (from the beginning of the simulation day) is used: top soil layer - infiltration is also counts in concentration calculation */

	for (layer = 0; layer < N_SOILLAYERS; layer++) soilwAVAIL[layer] = ws->soilw_pre[layer] - sprop->VWChw[layer] / sitec->soillayer_thickness[layer] / water_density;

	/*---------------------------------------------------------------------------------*/
	/* leaching fluxes for the 10 dissolving material types */
	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		/* special layers in case on GW: CFlayer and GWlayer */
		if (sprop->GWlayer != DATA_GAP && (layer == CFlayer || layer == GWlayer))
		{ 
			for (dm = 0; dm < N_DISSOLVMATER; dm++)
			{
				/* CFlayer*/
				if (layer == CFlayer && layer != GWlayer)
				{ 
					if (!errorCode && capillary_leaching(dm, sitec, sprop, soilInfo, ws, wf, &dismatLeachNORM, &dismatLeachCAPIL, &dischargeNORM, &dischargeCAPIL, &rechargeNORM, &rechargeCAPIL))
					{
						printf("ERROR in capillary_leaching.c from multilayer_leaching.c\n");
						errorCode = 1;
					}
				}	
				/* GWlayer*/
				else
				{	
					if (!errorCode && groundwater_leaching(dm, sitec, sprop, soilInfo, ws, wf, &dismatLeachNORM, &dismatLeachCAPIL, &dischargeNORM, &dischargeCAPIL, &rechargeNORM, &rechargeCAPIL))
					{
						printf("ERROR in groundwater_leaching.c from multilayer_leaching.c\n");
						errorCode = 1;
					}
				}

				soilInfo->dismatLeach[dm][layer] = dismatLeachNORM + dismatLeachCAPIL;

				soilInfo->dismatGWdischarge[dm][layer] = dischargeNORM + dischargeCAPIL;
				soilInfo->dismatGWrecharge[dm][layer] = rechargeNORM + rechargeCAPIL;
			}
		}
		else
		{
			wflux = wf->soilwFlux[layer];

			if (wflux != 0)
			{				
				for (dm = 0; dm < N_DISSOLVMATER; dm++)
				{
					conc0[dm] =  (soilInfo->contentDISSOLV_soil[dm][layer] / ws->soilw[layer]);
					if (layer < N_SOILLAYERS - 1)
					{		
						conc1[dm] =  (soilInfo->contentDISSOLV_soil[dm][layer + 1] / ws->soilw[layer + 1]);

						if (wflux > 0)
						{
							if (wflux > soilwAVAIL[layer]) wflux = soilwAVAIL[layer];
							soilInfo->dismatLeach[dm][layer] = conc0[dm] * wflux;
						}
						else
						{
							if (fabs(wflux) > soilwAVAIL[layer+1]) wflux = -1 * soilwAVAIL[layer+1];
							soilInfo->dismatLeach[dm][layer] = conc1[dm] * wflux;
						}
					}
					else
						soilInfo->dismatLeach[dm][layer] = conc0[dm] * wflux;
				}	

			}
		}
	}

	/* if GWD: soilwPercolDiffus from NORMcf layer */
	if (sprop->GWlayer != DATA_GAP && sprop->dz_NORMcf)
	{
		for (layer = CFlayer; layer < GWlayer; layer++)
		{ 
			for (dm = 0; dm < N_DISSOLVMATER; dm++)
			{
				soilInfo->dismatLeach_percolDiffus[dm][layer] = wf->soilwPercolDiffus_fromNORM[layer] * soilInfo->contentDISSOLV_NORMcf[dm] / sprop->soilw_NORMcf;
			}
		}
	}
	/*---------------------------------------------------------------------------------*/
	/* state update of pools: non-GWD (dismatLeach_fromAbove and dismatLeach) and GWD fluxes (percolDiffus, GWdischarge, GWdischarge) */

	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{	
		for (dm = 0; dm < N_DISSOLVMATER; dm++)
		{
			if (layer == 0)
				dismatLeach_fromAbove = 0;
			else
				dismatLeach_fromAbove = soilInfo->dismatLeach[dm][layer - 1];

	
			change = dismatLeach_fromAbove - soilInfo->dismatLeach[dm][layer] +
				soilInfo->dismatLeach_percolDiffus[dm][layer] + soilInfo->dismatGWdischarge[dm][layer] - soilInfo->dismatGWrecharge[dm][layer];

			soilInfo->contentDISSOLV_soil[dm][layer] += change;
	

			/* avoiding negative pool */
			if (soilInfo->contentDISSOLV_soil[dm][layer] < 0)
			{
				diff = soilInfo->contentDISSOLV_soil[dm][layer];
				if (soilInfo->dismatLeach[dm][layer] > fabs(diff))
					soilInfo->dismatLeach[dm][layer] += diff;
				else
				{
					if (soilInfo->dismatGWrecharge[dm][layer] > fabs(diff))
						soilInfo->dismatGWrecharge[dm][layer] += diff;
					else
					{
						printf("\n");
						printf("ERROR: negative content_CAPILgw in multilayer_leaching.c\n");
						errorCode = 1;
					}
				}

				
				/* correcting CFlayer-values if actual layer is CFlayer*/
				if (layer == CFlayer)
				{
					soilInfo->contentDISSOLV_NORMcf[dm] = 0;
					soilInfo->contentDISSOLV_CAPILcf[dm] = 0;
				}

				/* correcting CFlayer-values if bottom neighbourg is CFlayer*/
				if (layer + 1 == CFlayer && layer + 1 != GWlayer)
				{
					diffNORM = diff * (soilInfo->contentDISSOLV_NORMcf[dm] / (soilInfo->contentDISSOLV_NORMcf[dm] + soilInfo->contentDISSOLV_CAPILcf[dm]));
					diffCAPIL = diff - diffNORM;
					soilInfo->contentDISSOLV_NORMcf[dm] += diffNORM;
					soilInfo->contentDISSOLV_CAPILcf[dm] += diffCAPIL;
				}

				soilInfo->contentDISSOLV_soil[dm][layer] = 0;
			}
		}
	}


	/*---------------------------------------------------------------------------------*/
	/* calculation of src/snk variables   */
	for (dm = 0; dm < N_DISSOLVN; dm++)
	{		
		for (layer = 0; layer < N_SOILLAYERS; layer++)
		{
			ns->GWsnk_N += soilInfo->dismatGWrecharge[dm][layer];
			ns->GWsrc_N += soilInfo->dismatGWdischarge[dm][layer];
		}
	}
	for (dm = N_DISSOLVN; dm < N_DISSOLVMATER; dm++)
	{
		for (layer = 0; layer < N_SOILLAYERS; layer++)
		{
			cs->GWsnk_C += soilInfo->dismatGWrecharge[dm][layer];
			cs->GWsrc_C += soilInfo->dismatGWdischarge[dm][layer];
		}
	}

	/*---------------------------------------------------------------------------------*/
	/* deepleach calculation from the bottom layer */

	for (dm = 0; dm < N_DISSOLVN; dm++)              ns->Ndeepleach_snk += soilInfo->dismatLeach[dm][N_SOILLAYERS - 1];
	for (dm = N_DISSOLVN; dm < N_DISSOLVMATER; dm++) cs->Cdeepleach_snk += soilInfo->dismatLeach[dm][N_SOILLAYERS - 1];

	/*---------------------------------------------------------------------------------*/	
	/* transfer value: content_array ->NH4, NO3, DOC, DON  etc.*/

	/*  firsttime_flag=1 (after calculation note initial values, int partlyORtotal_flag=0 (BOUND or DISSOLV) is affected  */
	if (!errorCode && calc_DISSOLVandBOUND(1, 0, sprop, soilInfo))
	{
		printf("ERROR in calc_DISSOLVandBOUND.c for multilayer_leaching.c\n");
		errorCode = 1;
	}


	if (!errorCode && check_soilcontent(-1, 1, sprop, cs, ns, soilInfo))
	{
		printf("ERROR in check_soilcontent.c for multilayer_leaching.c\n");
		errorCode = 1;
	}

	
	
	return (errorCode);
}

