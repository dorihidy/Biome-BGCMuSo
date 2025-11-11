/* 
check_virtualLayer_balance.c
control of VWC values

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

int check_virtualLayer_balance(control_struct* ctrl, soilInfo_struct* soilInfo, soilprop_struct* sprop, wflux_struct* wf)
{


	int errorCode=0;
	int dm, GWlayer, CFlayer;
	double percolCF, percolGW;

	GWlayer = (int)sprop->GWlayer;
	CFlayer = (int)sprop->CFlayer;

	double change;

	/*--------------------------------------------------------------------------*/
	/* calculation of water balance in virtual layers */

	/* if CF layer is on top soil layers, no percolation from above layer, but evaporation counts if GW-layer is not in top soil layer */
	if (CFlayer == 0)
		percolCF = wf->infiltPOT;
	else
		percolCF = wf->soilwPercol[CFlayer - 1];

	if (GWlayer == 0)
		percolGW = wf->infiltPOT;
	else
		percolGW = wf->soilwPercol[GWlayer - 1];

	if (sprop->dz_NORMcf)
	{
		
		wf->inflow_NORMcf = percolCF + wf->soilwDiffus_aboveCFlayer_vs_NORMcf + wf->GWdischargeNORMcf + wf->FRZ_to_NORM[CFlayer];
		wf->inflow_CAPILcf = wf->soilwPercol_NORMvsCAPILcf + wf->soilwDiffus_aboveCFlayer_vs_CAPILcf + wf->GWdischargeCAPILcf + wf->soilwDiffus_NORMvsCAPILcf + wf->FRZ_to_CAPIL[CFlayer];
	}
	else
	{
		wf->inflow_NORMcf = 0;
		wf->inflow_CAPILcf = percolCF + wf->soilwPercol_NORMvsCAPILcf + wf->soilwDiffus_aboveCFlayer_vs_CAPILcf + wf->GWdischargeCAPILcf + wf->soilwDiffus_NORMvsCAPILcf + wf->FRZ_to_CAPIL[CFlayer];
	}

	if (sprop->dz_NORMgw)
	{
		wf->inflow_NORMgw = percolGW + wf->soilwDiffus_aboveGWlayer_vs_NORMgw + wf->GWdischargeNORMgw + wf->FRZ_to_NORM[GWlayer];
		wf->inflow_CAPILgw = wf->soilwPercol_NORMvsCAPILgw + wf->soilwDiffus_aboveGWlayer_vs_CAPILgw + wf->GWdischargeCAPILgw + wf->soilwDiffus_NORMvsCAPILgw + wf->soilwPercolDiffus_fromNORM[GWlayer] + wf->FRZ_to_CAPIL[GWlayer];
	}
	else
	{
		wf->inflow_NORMgw = 0;
		wf->inflow_CAPILgw = percolGW + wf->soilwPercol_NORMvsCAPILgw + wf->soilwDiffus_aboveGWlayer_vs_CAPILgw + wf->GWdischargeCAPILgw + wf->soilwDiffus_NORMvsCAPILgw + wf->soilwPercolDiffus_fromNORM[GWlayer] + +wf->FRZ_to_CAPIL[GWlayer];
	}


	wf->outflow_NORMcf = wf->soilwPercol_NORMvsCAPILcf + wf->soilwDiffusNORMcf + wf->soilwDiffus_NORMvsCAPILcf + wf->GWrecharge_NORMcf + wf->soilwPercolDiffus_fromNORM_total + wf->EVPsoilwNORMcf + wf->TRPsoilwNORMcf;
	wf->outflow_CAPILcf = wf->soilwPercol[CFlayer] + wf->soilwDiffusCAPILcf + wf->GWrecharge_CAPILcf + wf->EVPsoilwCAPILcf + wf->TRPsoilwCAPILcf;
	wf->outflow_NORMgw = wf->soilwPercol_NORMvsCAPILgw + wf->soilwDiffus_NORMvsCAPILgw + wf->GWrecharge_NORMgw + wf->TRPsoilw_NORMgw;
	wf->outflow_CAPILgw = wf->soilwPercol[GWlayer] + wf->GWrecharge_CAPILgw + wf->TRPsoilw_CAPILgw;


	if (ctrl->simyr == 0 && ctrl->yday == 0)
	{
		sprop->balance_NORMcf = 0;
		sprop->balance_CAPILcf = 0;
		sprop->balance_NORMgw = 0;
		sprop->balance_CAPILgw = 0;
	}
	else
	{
		if (sprop->dz_NORMcf) sprop->balance_NORMcf += sprop->soilw_NORMcf_pre + wf->inflow_NORMcf - wf->outflow_NORMcf - sprop->soilw_NORMcf;
		if (sprop->dz_CAPILcf) sprop->balance_CAPILcf += sprop->soilw_CAPILcf_pre + wf->inflow_CAPILcf - wf->outflow_CAPILcf - sprop->soilw_CAPILcf;
		if (sprop->dz_NORMgw) sprop->balance_NORMgw += sprop->soilw_NORMgw_pre + wf->inflow_NORMgw - wf->outflow_NORMgw - sprop->soilw_NORMgw;
		if (sprop->dz_CAPILgw) sprop->balance_CAPILgw += sprop->soilw_CAPILgw_pre + wf->inflow_CAPILgw - wf->outflow_CAPILgw - sprop->soilw_CAPILgw;

		if (fabs(sprop->balance_NORMcf) < CRIT_PREC_lenient && sprop->balance_NORMcf != 0) sprop->balance_NORMcf = 0;
		if (fabs(sprop->balance_CAPILcf) < CRIT_PREC_lenient && sprop->balance_CAPILcf != 0) sprop->balance_CAPILcf = 0;
		if (fabs(sprop->balance_NORMgw) < CRIT_PREC_lenient && sprop->balance_NORMgw != 0) sprop->balance_NORMgw = 0;
		if (fabs(sprop->balance_CAPILgw) < CRIT_PREC_lenient && sprop->balance_CAPILgw != 0) sprop->balance_CAPILgw = 0;
	}

	/*--------------------------------------------------------------------------*/
	

	for (dm = 0; dm < N_DISSOLVMATER; dm++)
	{
		change = soilInfo->dismatUNSATdecomp[dm][GWlayer] + soilInfo->dismatUNSATecofunc[dm][GWlayer] + soilInfo->dismatUNSATfertil[dm][GWlayer];

		if (GWlayer != 0) change += soilInfo->dismatLeach[dm][GWlayer - 1];
		change += soilInfo->dismatLeach_percolDiffus[dm][GWlayer];
		change += soilInfo->dismatGWdischarge[dm][GWlayer];
		change -= soilInfo->dismatLeach[dm][GWlayer];
		change -= soilInfo->dismatGWrecharge[dm][GWlayer];
	
		if (ctrl->simyr == 0 && ctrl->yday == 0)
			soilInfo->balance_UNSAT[dm] = 0;
		else
		{
			soilInfo->balance_UNSAT[dm] += soilInfo->content_zoneUNSATgw_pre[dm] + change - soilInfo->content_NORMgw[dm] - soilInfo->content_CAPILgw[dm];
			if (fabs(soilInfo->balance_UNSAT[dm]) < CRIT_PREC_lenient) soilInfo->balance_UNSAT[dm] = 0;

		}

		soilInfo->content_zoneUNSATgw_pre[dm] = soilInfo->content_NORMgw[dm] + soilInfo->content_CAPILgw[dm];

	}

	/*--------------------------------------------------------------------------*/

	for (dm = 0; dm < N_DISSOLVMATER; dm++)
	{
		if (sprop->soilw_NORMcf && soilInfo->content_NORMcf[dm] != 0 && soilInfo->content_NORMcf[dm] < CRIT_PREC_superRIG)
		{
			soilInfo->content_NORMcf[dm] = 0;
			if (dm >= N_DISSOLVinorgN && dm < N_DISSOLVN) if (soilInfo->content_NORMcf[dm + N_DISSOLVorgN] < CRIT_PREC_RIG) soilInfo->content_NORMcf[dm + N_DISSOLVorgN] = 0;
		}
		if (sprop->soilw_CAPILcf && soilInfo->content_CAPILcf[dm] != 0 && soilInfo->content_CAPILcf[dm] < CRIT_PREC_superRIG)
		{
			soilInfo->content_CAPILcf[dm] = 0;
			if (dm >= N_DISSOLVinorgN && dm < N_DISSOLVN) if (soilInfo->content_CAPILcf[dm + N_DISSOLVorgN] < CRIT_PREC_RIG) soilInfo->content_NORMgw[dm + N_DISSOLVorgN] = 0;
		}
		if (sprop->soilw_NORMgw && soilInfo->content_NORMgw[dm] != 0 && soilInfo->content_NORMgw[dm] < CRIT_PREC_superRIG)
		{
			soilInfo->content_NORMgw[dm] = 0;
			if (dm >= N_DISSOLVinorgN && dm < N_DISSOLVN) if (soilInfo->content_NORMgw[dm + N_DISSOLVorgN] < CRIT_PREC_RIG) soilInfo->content_NORMgw[dm + N_DISSOLVorgN] = 0;
		}
		if (sprop->soilw_CAPILgw && soilInfo->content_CAPILgw[dm] != 0 && soilInfo->content_CAPILgw[dm] < CRIT_PREC_superRIG)
		{
			soilInfo->content_CAPILgw[dm] = 0;
			if (dm >= N_DISSOLVinorgN && dm < N_DISSOLVN) if (soilInfo->content_CAPILgw[dm + N_DISSOLVorgN] < CRIT_PREC_RIG) soilInfo->content_CAPILgw[dm + N_DISSOLVorgN] = 0;
		}
	}

	return(errorCode);
}

