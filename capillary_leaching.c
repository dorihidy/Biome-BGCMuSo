/* 
groundwater_leaching.c
Calculating soil mineral nitrogen and DOC-DON leaching in capillary layers 
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

int capillary_leaching(int dm, const siteconst_struct* sitec,  soilprop_struct* sprop, soilInfo_struct* soilInfo, wstate_struct* ws, wflux_struct* wf,
	                   double* dismatLeachNORM, double* dismatLeachCAPIL, double* dischargeNORM, double* dischargeCAPIL, double* rechargeNORM, double* rechargeCAPIL)
{

	int errorCode, CFlayer, GWlayer, layer;

	double wflux,  diff, diff_total;
	double dismatLeach_NORMvsCAPIL, dismatLeach_NORMfromAbove, dismatLeach_CAPILfromAbove, percolDiffus_NORM, dismatLeach_fromAbove;
	double dismatLeachNORM_act, dismatLeachCAPIL_act, dischargeNORM_act, dischargeCAPIL_act, rechargeNORM_act, rechargeCAPIL_act;
	double soilw[N_SOILLAYERS], soilwAVAIL[N_SOILLAYERS];
	double soilwAVAIL_NORMcf, soilwAVAIL_CAPILcf, soilw_NORMcf, soilw_CAPILcf, soilw_CAPILgw, soilwAVAIL_CAPILgw;
	double change_NORM, change_CAPIL;

	errorCode = 0;

	dismatLeachNORM_act = wflux = dismatLeach_NORMvsCAPIL = dismatLeach_NORMfromAbove = dismatLeach_CAPILfromAbove = dischargeNORM_act = dischargeCAPIL_act = rechargeNORM_act = rechargeCAPIL_act = dismatLeachCAPIL_act = percolDiffus_NORM  = dismatLeach_fromAbove = 0;

	GWlayer = (int)sprop->GWlayer;
	CFlayer = (int)sprop->CFlayer;



	/*---------------------------------------------------------------------------------*/
	/* for concentration calculation original soilw data (from the beginning of the simulation day) is used: top soil layer - infiltration is also counts in concentration calculation */

	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		soilw[layer] = ws->soilw_pre[layer];
		soilwAVAIL[layer] = ws->soilw_pre[layer] - sprop->VWChw[layer] / sitec->soillayer_thickness[layer] / water_density;
	}
	soilw_NORMcf = sprop->soilw_NORMcf_pre;
	soilw_CAPILcf = sprop->soilw_CAPILcf_pre;
	soilw_CAPILgw = sprop->soilw_CAPILgw_pre;
	soilwAVAIL_NORMcf = sprop->soilw_NORMcf_pre - sprop->VWChw[CFlayer] / sprop->dz_NORMcf / water_density;
	soilwAVAIL_CAPILcf = sprop->soilw_CAPILcf_pre - sprop->VWChw[CFlayer] / sprop->dz_CAPILcf / water_density;
	soilwAVAIL_CAPILgw = sprop->soilw_CAPILgw_pre - sprop->VWChw[GWlayer] / sprop->dz_CAPILgw / water_density;
	/*---------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
	/*---------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
	/* flux between NORM and CAPIL */
	if (sprop->dz_NORMcf && sprop->dz_CAPILcf)
	{
		wflux = (wf->soilwPercol_NORMvsCAPILcf + wf->soilwDiffus_NORMvsCAPILcf);
		if (wflux > 0)
		{
			if (wflux > soilwAVAIL_NORMcf) wflux = soilwAVAIL_NORMcf;
			dismatLeach_NORMvsCAPIL = wflux * (soilInfo->contentDISSOLV_NORMcf[dm] / soilw_NORMcf);
		}
		else
		{
			if (fabs(wflux) > soilwAVAIL_CAPILcf) wflux = -1 * soilwAVAIL_CAPILcf;
			dismatLeach_NORMvsCAPIL = wflux * (soilInfo->contentDISSOLV_CAPILcf[dm] / soilw_CAPILcf);
		}
	}

	/*-------------------------------------------------------------------*/
	/* fluxes from above: top soil CFlayer is special */
	if (CFlayer != 0)
	{
		/* NORMcf */
		if (sprop->dz_NORMcf)
		{
			wflux = wf->soilwPercol[CFlayer - 1] + wf->soilwDiffus_aboveCFlayer_vs_NORMcf;
			if (wflux > 0)
			{
				if (wflux > soilwAVAIL[CFlayer - 1]) wflux = soilwAVAIL[CFlayer - 1];
				dismatLeach_NORMfromAbove = wflux * (soilInfo->contentDISSOLV_soil[dm][CFlayer - 1] / soilw[CFlayer - 1]);
			}
			else
			{
				if (fabs(wflux) > soilwAVAIL_NORMcf) wflux = -1 * soilwAVAIL_NORMcf;
				dismatLeach_NORMfromAbove = wflux * (soilInfo->contentDISSOLV_NORMcf[dm] / soilw_NORMcf);
			}
		}

		/* CAPILcf: if there is no NORMcf - percolation from upper CFlayer enter into CAPILcf first */
		if (sprop->dz_CAPILcf)
		{
			if (!sprop->dz_NORMcf)
				wflux = wf->soilwPercol[CFlayer - 1] + wf->soilwDiffus_aboveCFlayer_vs_CAPILcf;
			else
				wflux = wf->soilwDiffus_aboveCFlayer_vs_CAPILcf;

			if (wflux > 0)
			{
				if (wflux > soilwAVAIL[CFlayer - 1]) wflux = soilwAVAIL[CFlayer - 1];
				dismatLeach_CAPILfromAbove = wflux * (soilInfo->contentDISSOLV_soil[dm][CFlayer - 1] / soilw[CFlayer - 1]);
			}
			else
			{
				if (fabs(wflux) > soilwAVAIL_CAPILcf) wflux = soilwAVAIL_CAPILcf;
				dismatLeach_CAPILfromAbove = wflux * (soilInfo->contentDISSOLV_CAPILcf[dm] / soilw_CAPILcf);
			}
		}
		soilInfo->dismatLeach[dm][CFlayer - 1] = dismatLeach_NORMfromAbove + dismatLeach_CAPILfromAbove;

	}
	else
	{
		dismatLeach_NORMfromAbove = 0;
		dismatLeach_CAPILfromAbove = 0;

	}

	/*-------------------------------------------------------------------*/
	/* leaching fluxes (percol+diffus) - mixed CFlayer! */

	/* NORMcf */
	if (sprop->dz_NORMcf)
	{
		/* if there is not CAPILcf - percolation from CF layer is from NORMcf*/
		if (sprop->dz_CAPILcf)
			wflux = wf->soilwDiffusNORMcf;
		else
			wflux = wf->soilwDiffusNORMcf + wf->soilwPercol[CFlayer];

		if (wflux > 0)
		{
			if (wflux > soilwAVAIL_NORMcf) wflux = soilwAVAIL_NORMcf;
			dismatLeachNORM_act = wflux * (soilInfo->contentDISSOLV_NORMcf[dm] / soilw_NORMcf);
		}
		else
		{
			if (CFlayer + 1 != GWlayer || sprop->dz_CAPILgw == 0)
			{
				if (fabs(wflux) > soilwAVAIL[CFlayer + 1]) wflux = -1 * soilwAVAIL[CFlayer + 1];
				dismatLeachNORM_act = wflux * (soilInfo->contentDISSOLV_soil[dm][CFlayer + 1] / ws->soilw[CFlayer + 1]);
			}
			else
			{
				if (fabs(wflux) > soilwAVAIL_CAPILgw) wflux = -1 * soilwAVAIL_CAPILgw;
				dismatLeachNORM_act = wflux * (soilInfo->contentDISSOLV_CAPILgw[dm] / soilw_CAPILgw);
			}
		}
	}

	/* CAPILcf */
	if (sprop->dz_CAPILcf)
	{ 
		wflux = wf->soilwPercol[CFlayer] + wf->soilwDiffusCAPILcf;

		if (wflux >= 0)
		{
			if (wflux > soilwAVAIL_CAPILcf) wflux = soilwAVAIL_CAPILcf;
			dismatLeachCAPIL_act = wflux * (soilInfo->contentDISSOLV_CAPILcf[dm] / soilw_CAPILcf);
		}
		else
		{
			/* if bottom neigbourh is GWlayer - CAPILgw is the other layer*/
			if (CFlayer + 1 != GWlayer || sprop->dz_CAPILgw == 0)
			{
				if (fabs(wflux) > ws->soilw[CFlayer + 1]) wflux = -1 * ws->soilw[CFlayer + 1];
				dismatLeachCAPIL_act = wflux * (soilInfo->contentDISSOLV_soil[dm][CFlayer + 1] / ws->soilw[CFlayer + 1]);
			}
			else
			{
				if (fabs(wflux) > soilwAVAIL_CAPILgw) wflux = -1 * soilwAVAIL_CAPILgw;
				dismatLeachCAPIL_act = wflux * (soilInfo->contentDISSOLV_CAPILgw[dm] / soilw_CAPILgw);
			}
		}
	}
	/* GWdischarge  */
	dischargeNORM_act = wf->GWdischargeNORMcf * soilInfo->GWconc[dm];
	dischargeCAPIL_act = wf->GWdischargeCAPILcf * soilInfo->GWconc[dm];

	/* GWrecharge: limitation of flux - size of the pool */
	if (sprop->dz_NORMcf)
	{
		wflux = wf->GWrecharge_NORMcf;
		if (wflux > soilwAVAIL_NORMcf) wflux = soilwAVAIL_NORMcf;
		rechargeNORM_act = wflux * (soilInfo->contentDISSOLV_NORMcf[dm] / soilw_NORMcf);
	}
	else
		rechargeNORM_act = 0;

	if (sprop->dz_CAPILcf)
	{
		wflux = wf->GWrecharge_CAPILcf;
		if (wflux > soilwAVAIL_CAPILcf) wflux = soilwAVAIL_CAPILcf;
		rechargeCAPIL_act = wflux * (soilInfo->contentDISSOLV_CAPILcf[dm] / soilw_CAPILcf);
	}
	else
		rechargeCAPIL_act = 0;

	/* PercolDiffus */
	if (sprop->dz_NORMcf)
	{
		wflux = wf->soilwPercolDiffus_fromNORM_total;
		if (wflux > soilwAVAIL_NORMcf) wflux = soilwAVAIL_NORMcf;
		percolDiffus_NORM = wf->soilwPercolDiffus_fromNORM_total * (soilInfo->contentDISSOLV_NORMcf[dm] / soilw_NORMcf);
	}

	/* state update */
	change_NORM = dismatLeach_NORMfromAbove + dischargeNORM_act - dismatLeach_NORMvsCAPIL - percolDiffus_NORM - rechargeNORM_act - dismatLeachNORM_act;
	change_CAPIL = dismatLeach_CAPILfromAbove + dischargeCAPIL_act + dismatLeach_NORMvsCAPIL + percolDiffus_NORM - rechargeCAPIL_act - dismatLeachCAPIL_act;

	soilInfo->contentDISSOLV_NORMcf[dm] += change_NORM;
	soilInfo->contentDISSOLV_CAPILcf[dm] += change_CAPIL;

	if (soilInfo->contentDISSOLV_NORMcf[dm] < 0 && fabs(soilInfo->contentDISSOLV_NORMcf[dm]) < CRIT_PREC) soilInfo->contentDISSOLV_NORMcf[dm] = 0;
	if (soilInfo->contentDISSOLV_CAPILcf[dm] < 0 && fabs(soilInfo->contentDISSOLV_CAPILcf[dm]) < CRIT_PREC) soilInfo->contentDISSOLV_CAPILcf[dm] = 0;

	/* avoiding negative pool: propoertionally reduce the sink flows */
	if (soilInfo->contentDISSOLV_NORMcf[dm] < 0)
	{
		diff = soilInfo->contentDISSOLV_NORMcf[dm];
		soilInfo->contentDISSOLV_NORMcf[dm] -= change_NORM;
		soilInfo->contentDISSOLV_CAPILcf[dm] -= change_CAPIL;

		diff_total = rechargeNORM_act;

		if (dismatLeachNORM_act > 0) diff_total += dismatLeachNORM_act;
		if (dismatLeach_NORMvsCAPIL > 0) diff_total += dismatLeach_NORMvsCAPIL;
		if (dismatLeach_NORMfromAbove < 0) diff_total += fabs(dismatLeach_NORMfromAbove);


		rechargeNORM_act += diff * rechargeNORM_act / diff_total;
		if (dismatLeachNORM_act > 0) dismatLeachNORM_act += diff * dismatLeachNORM_act / diff_total;
		if (dismatLeach_NORMvsCAPIL > 0) dismatLeach_NORMvsCAPIL += diff * dismatLeach_NORMvsCAPIL / diff_total;
		if (dismatLeach_NORMfromAbove < 0)
		{
			dismatLeach_NORMfromAbove += diff * dismatLeach_NORMfromAbove / diff_total;
			soilInfo->dismatLeach[dm][CFlayer - 1] += diff * dismatLeach_NORMfromAbove / diff_total;
			soilInfo->contentDISSOLV_soil[dm][CFlayer - 1] -= diff * dismatLeach_NORMfromAbove / diff_total;
		}

		soilInfo->dismatLeach[dm][CFlayer - 1] = dismatLeach_NORMfromAbove + dismatLeach_CAPILfromAbove;

		change_NORM = dismatLeach_NORMfromAbove + dischargeNORM_act - dismatLeach_NORMvsCAPIL - percolDiffus_NORM - rechargeNORM_act - dismatLeachNORM_act;
		change_CAPIL = dismatLeach_CAPILfromAbove + dischargeCAPIL_act + dismatLeach_NORMvsCAPIL + percolDiffus_NORM - rechargeCAPIL_act - dismatLeachCAPIL_act;

		soilInfo->contentDISSOLV_NORMcf[dm] += change_NORM;
		soilInfo->contentDISSOLV_CAPILcf[dm] += change_CAPIL;

		if (soilInfo->contentDISSOLV_NORMcf[dm] < 0 && fabs(soilInfo->contentDISSOLV_NORMcf[dm]) < CRIT_PREC) soilInfo->contentDISSOLV_NORMcf[dm] = 0;
		if (soilInfo->contentDISSOLV_CAPILcf[dm] < 0 && fabs(soilInfo->contentDISSOLV_CAPILcf[dm]) < CRIT_PREC) soilInfo->contentDISSOLV_CAPILcf[dm] = 0;

		if (soilInfo->contentDISSOLV_NORMcf[dm] < 0)
		{
			printf("\n");
			printf("ERROR: negative contentDISSOLV_CAPILcf in multilayer_leaching.c\n");
			errorCode = 1;
		}

	}

	/* avoiding negative pool: propoertionally reduce the sink flows */
	if (soilInfo->contentDISSOLV_CAPILcf[dm] < 0)
	{
		diff = soilInfo->contentDISSOLV_CAPILcf[dm];
		soilInfo->contentDISSOLV_NORMcf[dm] -= change_NORM;
		soilInfo->contentDISSOLV_CAPILcf[dm] -= change_CAPIL;

		diff_total = rechargeCAPIL_act;

		if (dismatLeachCAPIL_act > 0) diff_total += dismatLeachCAPIL_act;
		if (dismatLeach_NORMvsCAPIL < 0) diff_total += fabs(dismatLeach_NORMvsCAPIL);
		if (dismatLeach_CAPILfromAbove < 0) diff_total += fabs(dismatLeach_CAPILfromAbove);


		rechargeCAPIL_act += diff * rechargeCAPIL_act / diff_total;
		if (dismatLeachCAPIL_act > 0) dismatLeachCAPIL_act += diff * dismatLeachCAPIL_act / diff_total;
		if (dismatLeach_NORMvsCAPIL < 0) dismatLeach_NORMvsCAPIL += diff * dismatLeach_NORMvsCAPIL / diff_total;
		if (dismatLeach_CAPILfromAbove < 0)
		{
			dismatLeach_CAPILfromAbove += diff * dismatLeach_CAPILfromAbove / diff_total;
			soilInfo->dismatLeach[dm][CFlayer - 1] += diff * dismatLeach_CAPILfromAbove / diff_total;
			soilInfo->contentDISSOLV_soil[dm][CFlayer - 1] -= diff * dismatLeach_CAPILfromAbove / diff_total;
		}

		soilInfo->dismatLeach[dm][CFlayer - 1] = dismatLeach_NORMfromAbove + dismatLeach_CAPILfromAbove;

		change_NORM = dismatLeach_NORMfromAbove + dischargeNORM_act - dismatLeach_NORMvsCAPIL - percolDiffus_NORM - rechargeNORM_act - dismatLeachNORM_act;
		change_CAPIL = dismatLeach_CAPILfromAbove + dischargeCAPIL_act + dismatLeach_NORMvsCAPIL + percolDiffus_NORM - rechargeCAPIL_act - dismatLeachCAPIL_act;

		soilInfo->contentDISSOLV_NORMcf[dm] += change_NORM;
		soilInfo->contentDISSOLV_CAPILcf[dm] += change_CAPIL;

		if (soilInfo->contentDISSOLV_NORMcf[dm] < 0 && fabs(soilInfo->contentDISSOLV_NORMcf[dm]) < CRIT_PREC) soilInfo->contentDISSOLV_NORMcf[dm] = 0;
		if (soilInfo->contentDISSOLV_CAPILcf[dm] < 0 && fabs(soilInfo->contentDISSOLV_CAPILcf[dm]) < CRIT_PREC) soilInfo->contentDISSOLV_CAPILcf[dm] = 0;

		if (soilInfo->contentDISSOLV_CAPILcf[dm] < 0)
		{
			printf("\n");
			printf("ERROR: negative contentDISSOLV_CAPILcf in multilayer_leaching.c\n");
			errorCode = 1;
		}

	}

	*dismatLeachNORM = dismatLeachNORM_act;
	*dismatLeachCAPIL = dismatLeachCAPIL_act;
	*dischargeNORM = dischargeNORM_act;
	*dischargeCAPIL = dischargeCAPIL_act;
	*rechargeNORM = rechargeNORM_act;
	*rechargeCAPIL = rechargeCAPIL_act;

	return (errorCode);
}

