/*
state_update.c
Resolve the fluxes in bgc.c daily loop to update state variables

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGCMuSo v7.0.
Original code: Copyright 2000, Peter E. Thornton
Numerical Terradynamic Simulation Group, The University of Montana, USA
Modified code: Copyright 2025, D. Hidy [dori.hidy@gmail.com]
Hungarian Academy of Sciences, Hungary
See the website of Biome-BGCMuSo at http://nimbus.elte.hu/bbgc/ for documentation, model executable and example input files.
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Modified:
4/17/2000 (PET): Included the new denitrification flux. See daily_allocation.c
for complete description of this change.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include "ini.h"
#include "bgc_struct.h"
#include "bgc_func.h"
#include "bgc_constants.h"

int water_state_update(wflux_struct* wf, wstate_struct* ws)
{
	/* daily update of the water state variables */

	int errorCode = 0;
	int layer;


	/* snoww */
	ws->snoww += wf->prcp_to_snoww;
	ws->snoww -= wf->snoww_to_soilw;
	ws->snoww -= wf->SUBLsnoww;

	/* precipitation fluxes */
	ws->canopyw += wf->prcp_to_canopyw;
	ws->canopyw -= wf->EVPcanopyw;
	ws->canopyw -= wf->canopyw_to_soilw;


	/* precipitation src */
	ws->prcp_src += wf->prcp_to_canopyw;
	ws->prcp_src += (wf->prcp_to_soilSurface - wf->IRG_to_prcp);
	ws->prcp_src += wf->prcp_to_snoww;


	/* evapotranspiration/sublimation snk */
	ws->EVPcanopyw_snk += wf->EVPcanopyw;
	ws->snowSUBL_snk += wf->SUBLsnoww;
	ws->soilEVP_snk += wf->EVPsoilw;
	ws->TRP_snk += wf->TRPsoilw_SUM;
	ws->pondEVP_snk += wf->EVPpondw;


	/* runoff */
	ws->runoff_snk += wf->pondw_to_runoff;

	ws->deeppercolation_snk += wf->soilwFlux[N_SOILLAYERS - 1];

	/* groundwater src/snk */
	wf->GWrecharge_total = wf->GWrecharge_CAPILcf + wf->GWrecharge_CAPILgw + wf->GWrecharge_NORMcf + wf->GWrecharge_NORMgw + wf->GWrecharge_lastCAPIL;
	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		wf->TRPfromGW_total    += wf->TRPfromGW[layer];
		wf->GWdischarge_total += wf->GWdischarge[layer];
		wf->soilwPercolDiffus_fromNORM_total += wf->soilwPercolDiffus_fromNORM[layer];
	}

	ws->GWsrc_W += wf->EVPfromGW;
	ws->GWsrc_W += wf->GW_to_pondw;
	ws->GWsrc_W += wf->TRPfromGW_total;
	ws->GWsrc_W += wf->GWdischarge_total;
	ws->GWsnk_W += wf->GWrecharge_total;


	/* flooding src */
	ws->FLsrc_W += wf->FL_to_soilwTOTAL;

	/* irrigating src*/
	ws->IRGsrc_W += wf->IRG_to_prcp;
	ws->IRGsrc_W += wf->IRG_to_soilw;
	ws->IRGsrc_W += wf->IRG_to_soilSurface;
	
	return (errorCode);
}

int CN_state_update(const siteconst_struct* sitec, const epconst_struct* epc, soilInfo_struct* soilInfo, soilprop_struct* sprop, control_struct* ctrl, epvar_struct* epv,
	                cflux_struct* cf, nflux_struct* nf, cstate_struct* cs, nstate_struct* ns, int alloc, int evergreen)
{
	/* daily update of the carbon state variables */

	int errorCode = 0;
	int layer, pp, dm, GWlayer, CFlayer;
	double leafc_to_litr, leafn_to_litr, frootc_to_litr, frootn_to_litr, yieldc_to_litr, yieldn_to_litr, softstemc_to_litr, softstemn_to_litr;
	double ratioNORM, ratioCAPIL, ratioSAT, diffNORM, diffCAPIL, diffSAT;

	/* C state variables are updated below in the order of the relevant fluxes in the daily model loop */

	/* NOTE: Mortality fluxes are all accounted for in a separate routine, which is to be called after this routine.
	This is a special case where the updating of state variables is order-sensitive, since otherwise the complications of possibly having
	mortality fluxes drive the pools negative would create big, unnecessary headaches. */

	/* 0. Initialization of local variables */

	GWlayer = (int)sprop->GWlayer;
	CFlayer = (int)sprop->CFlayer;

	/* *****************************************************************************************************/
	/* 0. firsttime_flag=0 (before calculation note initial values, partlyORtotal_flag=1 (TOTAL (BOUND+DISSOLV) is affected  */

	if (!errorCode && calc_DISSOLVandBOUND(0, 1, sprop, soilInfo))
	{
		printf("ERROR in calc_DISSOLVandBOUND.c for state_update.c\n");
		errorCode = 1;
	}

	/***************************************************************************************************************************************************/
	/* 1. Phenology fluxes */
	if (!errorCode && epc->leaf_cn && CNratio_control(cs, epc->leaf_cn, cs->leafc, ns->leafn, cf->leafc_transfer_to_leafc, nf->leafn_transfer_to_leafn, 0)) 
	{
		cs->leafc               += cf->leafc_transfer_to_leafc;
		cs->leafc_transfer      -= cf->leafc_transfer_to_leafc;
		ns->leafn               = cs->leafc / epc->leaf_cn;
		ns->leafn_transfer      = cs->leafc_transfer / epc->leaf_cn;
	}
	else if (epc->leaf_cn)		
	{
		if (!errorCode) printf("ERROR in leaf CN calculation in state_update.c\n");
		errorCode=1;
	}

	if (!errorCode && epc->froot_cn && CNratio_control(cs, epc->froot_cn, cs->frootc, ns->frootn, cf->frootc_transfer_to_frootc, nf->frootn_transfer_to_frootn, 0)) 
	{
		cs->frootc               += cf->frootc_transfer_to_frootc;
		cs->frootc_transfer      -= cf->frootc_transfer_to_frootc;
		ns->frootn               = cs->frootc / epc->froot_cn;
		ns->frootn_transfer      = cs->frootc_transfer / epc->froot_cn;
	}
	else if (epc->froot_cn)		
	{
		if (!errorCode) printf("ERROR in froot CN calculation in state_update.c\n");
		errorCode=1;
	}

	if (!errorCode && epc->yield_cn && CNratio_control(cs, epc->yield_cn, cs->yieldc, ns->yieldn, cf->yieldc_transfer_to_yield, nf->yieldn_transfer_to_yieldn, 0)) 
	{
		cs->yieldc               += cf->yieldc_transfer_to_yield;
		cs->yieldc_transfer      -= cf->yieldc_transfer_to_yield;
		ns->yieldn               = cs->yieldc / epc->yield_cn;
		ns->yieldn_transfer      = cs->yieldc_transfer / epc->yield_cn;
	}
	else if (epc->yield_cn)	
	{	
		if (!errorCode) printf("ERROR in yield CN calculation in state_update.c\n");
		errorCode=1;
	}

	if (!errorCode && epc->softstem_cn && CNratio_control(cs, epc->softstem_cn, cs->softstemc, ns->softstemn, cf->softstemc_transfer_to_softstemc, nf->softstemn_transfer_to_softstemn, 0)) 
	{
		cs->softstemc               += cf->softstemc_transfer_to_softstemc;
		cs->softstemc_transfer      -= cf->softstemc_transfer_to_softstemc;
		ns->softstemn               = cs->softstemc / epc->softstem_cn;
		ns->softstemn_transfer      = cs->softstemc_transfer / epc->softstem_cn;
	}
	else if (epc->softstem_cn)		
	{
		if (!errorCode) printf("ERROR in softstem CN calculation in state_update.c\n");
		errorCode=1;
	}
	
	if (!errorCode && epc->livewood_cn && CNratio_control(cs, epc->livewood_cn, cs->livestemc, ns->livestemn, cf->livestemc_transfer_to_livestemc, nf->livestemn_transfer_to_livestemn, 0)) 
	{
		cs->livestemc               += cf->livestemc_transfer_to_livestemc;
		cs->livestemc_transfer      -= cf->livestemc_transfer_to_livestemc;
		ns->livestemn               = cs->livestemc / epc->livewood_cn;
		ns->livestemn_transfer      = cs->livestemc_transfer / epc->livewood_cn;
	}
	else if (epc->livewood_cn)		
	{
		if (!errorCode) printf("ERROR in livestem CN calculation in state_update.c\n");
		errorCode=1;
	}

	if (!errorCode && epc->livewood_cn && CNratio_control(cs, epc->livewood_cn, cs->livecrootc, ns->livecrootn, cf->livecrootc_transfer_to_livecrootc, nf->livecrootn_transfer_to_livecrootn, 0)) 
	{
		cs->livecrootc               += cf->livecrootc_transfer_to_livecrootc;
		cs->livecrootc_transfer      -= cf->livecrootc_transfer_to_livecrootc;
		ns->livecrootn               = cs->livecrootc / epc->livewood_cn;
		ns->livecrootn_transfer      = cs->livecrootc_transfer / epc->livewood_cn;
	}
	else if (epc->livewood_cn)		
	{
		if (!errorCode) printf("ERROR in livecroot CN calculation in state_update.c\n");
		errorCode=1;
	}

	if (!errorCode && epc->deadwood_cn && CNratio_control(cs, epc->deadwood_cn, cs->deadstemc, ns->deadstemn, cf->deadstemc_transfer_to_deadstemc, nf->deadstemn_transfer_to_deadstemn, 0)) 
	{
		cs->deadstemc               += cf->deadstemc_transfer_to_deadstemc;
		cs->deadstemc_transfer      -= cf->deadstemc_transfer_to_deadstemc;
		ns->deadstemn               = cs->deadstemc / epc->deadwood_cn;
		ns->deadstemn_transfer      = cs->deadstemc_transfer / epc->deadwood_cn;
	}
	else if (epc->deadwood_cn)		
	{
		if (!errorCode) printf("ERROR in deadstem CN calculation in state_update.c\n");
		errorCode=1;
	}

	if (!errorCode && epc->deadwood_cn && CNratio_control(cs, epc->deadwood_cn, cs->deadcrootc, ns->deadcrootn, cf->deadcrootc_transfer_to_deadcrootc, nf->deadcrootn_transfer_to_deadcrootn, 0)) 
	{
		cs->deadcrootc               += cf->deadcrootc_transfer_to_deadcrootc;
		cs->deadcrootc_transfer      -= cf->deadcrootc_transfer_to_deadcrootc;
		ns->deadcrootn               = cs->deadcrootc / epc->deadwood_cn;
		ns->deadcrootn_transfer      = cs->deadcrootc_transfer / epc->deadwood_cn;
	}
	else if (epc->deadwood_cn)		
	{
		if (!errorCode) printf("ERROR in deadcroot CN calculation in state_update.c\n");
		errorCode=1;
	}

	/***************************************************************************************************************************************************/
	/* 2. Aboveground pool litterfall and retranslocation to the first soil layer */
	leafc_to_litr = cf->leafc_to_litr1c + cf->leafc_to_litr2c + cf->leafc_to_litr3c + cf->leafc_to_litr4c;
	leafn_to_litr = nf->leafn_to_litr1n + nf->leafn_to_litr2n + nf->leafn_to_litr3n + nf->leafn_to_litr4n;
		
	/* litter turns into the first three soil layers  (non-woody biomass: proportion to soil layer thickness, woody-biomass: higher propotion in layer2 */
	
	if (!errorCode && epc->leaf_cn && CNratio_control(cs, epc->leaf_cn, cs->leafc, ns->leafn, leafc_to_litr, leafn_to_litr, epc->leaflitr_cn)) 
	{
		for (layer = 0; layer < N_SOILLAYERS; layer++)
		{
			cs->litr1c[layer] += (cf->leafc_to_litr1c) * sprop->PROPlayerDC[layer];
			cs->litr2c[layer] += (cf->leafc_to_litr2c) * sprop->PROPlayerDC[layer];
			cs->litr3c[layer] += (cf->leafc_to_litr3c) * sprop->PROPlayerDC[layer];
			cs->litr4c[layer] += (cf->leafc_to_litr4c) * sprop->PROPlayerDC[layer];

			ns->litr1n[layer] += (nf->leafn_to_litr1n) * sprop->PROPlayerDC[layer];
			ns->litr2n[layer] += (nf->leafn_to_litr2n) * sprop->PROPlayerDC[layer];
			ns->litr3n[layer] += (nf->leafn_to_litr3n) * sprop->PROPlayerDC[layer];
			ns->litr4n[layer] += (nf->leafn_to_litr4n) * sprop->PROPlayerDC[layer];
		}

	
		cs->leafc      -= (cf->leafc_to_litr1c + cf->leafc_to_litr2c + cf->leafc_to_litr3c + cf->leafc_to_litr4c);
		ns->leafn       = cs->leafc / epc->leaf_cn;

		ns->retransn   += nf->leafn_to_retransn;   

		
	}
	else if (epc->leaflitr_cn || epc->leaf_cn)		
	{
		if (!errorCode) printf("ERROR in leaf_to_litr CN calculation in state_update.c\n");
		errorCode=1;
	}
	
	yieldc_to_litr = cf->yieldc_to_litr1c + cf->yieldc_to_litr2c + cf->yieldc_to_litr3c + cf->yieldc_to_litr4c;
	yieldn_to_litr = nf->yieldn_to_litr1n + nf->yieldn_to_litr2n + nf->yieldn_to_litr3n + nf->yieldn_to_litr4n;
	if (!errorCode && epc->yield_cn && CNratio_control(cs, epc->yield_cn, cs->yieldc, ns->yieldn, yieldc_to_litr, yieldn_to_litr, 0)) 
	{
		for (layer = 0; layer < N_SOILLAYERS; layer++)
		{
			cs->litr1c[layer] += (cf->yieldc_to_litr1c) * sprop->PROPlayerDC[layer];
			cs->litr2c[layer] += (cf->yieldc_to_litr2c) * sprop->PROPlayerDC[layer];
			cs->litr3c[layer] += (cf->yieldc_to_litr3c) * sprop->PROPlayerDC[layer];
			cs->litr4c[layer] += (cf->yieldc_to_litr4c) * sprop->PROPlayerDC[layer];

			ns->litr1n[layer] += (nf->yieldn_to_litr1n) * sprop->PROPlayerDC[layer];
			ns->litr2n[layer] += (nf->yieldn_to_litr2n) * sprop->PROPlayerDC[layer];
			ns->litr3n[layer] += (nf->yieldn_to_litr3n) * sprop->PROPlayerDC[layer];
			ns->litr4n[layer] += (nf->yieldn_to_litr4n) * sprop->PROPlayerDC[layer];
		}

	
		cs->yieldc -= (cf->yieldc_to_litr1c + cf->yieldc_to_litr2c + cf->yieldc_to_litr3c + cf->yieldc_to_litr4c);
		ns->yieldn      = cs->yieldc / epc->yield_cn;
	}
	else if (epc->yield_cn)		
	{
		if (!errorCode) printf("ERROR in yieldc_to_litr CN calculation in state_update.c\n");
		errorCode=1;
	}

	softstemc_to_litr = cf->softstemc_to_litr1c + cf->softstemc_to_litr2c + cf->softstemc_to_litr3c + cf->softstemc_to_litr4c;
	softstemn_to_litr = nf->softstemn_to_litr1n + nf->softstemn_to_litr2n + nf->softstemn_to_litr3n + nf->softstemn_to_litr4n;
	if (!errorCode && epc->softstem_cn && CNratio_control(cs, epc->softstem_cn, cs->softstemc, ns->softstemn, softstemc_to_litr, softstemn_to_litr, 0)) 
	{
		for (layer = 0; layer < N_SOILLAYERS; layer++)
		{

			cs->litr1c[layer] += (cf->softstemc_to_litr1c) * sprop->PROPlayerDC[layer];
			cs->litr2c[layer] += (cf->softstemc_to_litr2c) * sprop->PROPlayerDC[layer];
			cs->litr3c[layer] += (cf->softstemc_to_litr3c) * sprop->PROPlayerDC[layer];
			cs->litr4c[layer] += (cf->softstemc_to_litr4c) * sprop->PROPlayerDC[layer];

			ns->litr1n[layer] += (nf->softstemn_to_litr1n) * sprop->PROPlayerDC[layer];
			ns->litr2n[layer] += (nf->softstemn_to_litr2n) * sprop->PROPlayerDC[layer];
			ns->litr3n[layer] += (nf->softstemn_to_litr3n) * sprop->PROPlayerDC[layer];
			ns->litr4n[layer] += (nf->softstemn_to_litr4n) * sprop->PROPlayerDC[layer];
		}
				
	
		cs->softstemc  -= (cf->softstemc_to_litr1c + cf->softstemc_to_litr2c + cf->softstemc_to_litr3c + cf->softstemc_to_litr4c);

		ns->softstemn   = cs->softstemc / epc->softstem_cn;
	}
	else if (epc->softstem_cn)		
	{
		if (!errorCode) printf("ERROR in softstem_to_litr CN calculation in state_update.c\n");
		errorCode=1;
	}

	/***************************************************************************************************************************************************/
	/* 3. Belowground litterfall is distributed between the different soil layers */
	
	frootc_to_litr = cf->frootc_to_litr1c + cf->frootc_to_litr2c + cf->frootc_to_litr3c + cf->frootc_to_litr4c;
	frootn_to_litr = nf->frootn_to_litr1n + nf->frootn_to_litr2n + nf->frootn_to_litr3n + nf->frootn_to_litr4n;
	if (!errorCode && epc->froot_cn && CNratio_control(cs, epc->froot_cn, cs->frootc, ns->frootn, frootc_to_litr, frootn_to_litr, 0)) 
	{
		for (layer = 0; layer < N_SOILLAYERS; layer++)
		{
			cs->litr1c[layer] += cf->frootc_to_litr1c * epv->rootlengthProp[layer];
			cs->litr2c[layer] += cf->frootc_to_litr2c * epv->rootlengthProp[layer];
			cs->litr3c[layer] += cf->frootc_to_litr3c * epv->rootlengthProp[layer];
			cs->litr4c[layer] += cf->frootc_to_litr4c * epv->rootlengthProp[layer];

			ns->litr1n[layer] += nf->frootn_to_litr1n * epv->rootlengthProp[layer];
			ns->litr2n[layer] += nf->frootn_to_litr2n * epv->rootlengthProp[layer];
			ns->litr3n[layer] += nf->frootn_to_litr3n * epv->rootlengthProp[layer];
			ns->litr4n[layer] += nf->frootn_to_litr4n * epv->rootlengthProp[layer];
		}
		cs->frootc       -= frootc_to_litr;
		ns->frootn	      = cs->frootc / epc->froot_cn;
	}
	else if (epc->froot_cn)		
	{
		if (!errorCode) printf("ERROR in froot_to_litr CN calculation in state_update.c\n");
		errorCode=1;
	}

	/***************************************************************************************************************************************************/
	 /* 4. Livewood turnover fluxes */
	cs->deadstemc  += cf->livestemc_to_deadstemc;
	cs->livestemc  -= cf->livestemc_to_deadstemc;
	cs->deadcrootc += cf->livecrootc_to_deadcrootc;
	cs->livecrootc -= cf->livecrootc_to_deadcrootc;
	ns->retransn      += nf->livestemn_to_retransn;  
	ns->livestemn     -= nf->livestemn_to_retransn;
	ns->retransn      += nf->livecrootn_to_retransn; 
	ns->livecrootn    -= nf->livecrootn_to_retransn;
	

	ns->deadstemn  += nf->livestemn_to_deadstemn;
	ns->livestemn  -= nf->livestemn_to_deadstemn;
	ns->deadcrootn += nf->livecrootn_to_deadcrootn;
	ns->livecrootn -= nf->livecrootn_to_deadcrootn;
	
	/***************************************************************************************************************************************************/
	/* 6. Photosynthesis fluxes */
	cs->psnsun_src   += cf->psnsun_to_cpool;
	cs->psnshade_src += cf->psnshade_to_cpool;
	cs->cpool        += (cf->psnsun_to_cpool + cf->psnshade_to_cpool);

	/***************************************************************************************************************************************************/
	/* 7. Plant allocation flux, from N retrans pool */
	ns->npool		    += nf->retransn_to_npool_total;
	ns->retransn        -= nf->retransn_to_npool_total;

	/***************************************************************************************************************************************************/
	/* 8. Litter decomposition fluxes - MULTILAYER SOIL */
	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		/* Fluxes out of coarse woody debris into litter pools */
		cs->litr2c[layer] += cf->cwdc_to_litr2c[layer];
		cs->litr3c[layer] += cf->cwdc_to_litr3c[layer];
		cs->litr4c[layer] += cf->cwdc_to_litr4c[layer];
		cs->cwdc[layer] -= (cf->cwdc_to_litr2c[layer] + cf->cwdc_to_litr3c[layer] + cf->cwdc_to_litr4c[layer]);
		/* Fluxes out of labile litter pool */
		cs->HRlitr1_snk += cf->litr1_hr[layer];
		cs->litr1c[layer] -= cf->litr1_hr[layer];
		cs->soil1c[layer] += cf->litr1c_to_soil1c[layer];
		cs->litr1c[layer] -= cf->litr1c_to_soil1c[layer];
		/* Fluxes out of cellulose litter pool */
		cs->HRlitr2_snk += cf->litr2_hr[layer];
		cs->litr2c[layer] -= cf->litr2_hr[layer];
		cs->soil2c[layer] += cf->litr2c_to_soil2c[layer];
		cs->litr2c[layer] -= cf->litr2c_to_soil2c[layer];
		/* Fluxes from shielded to unshielded cellulose pools */
		cs->litr2c[layer] += cf->litr3c_to_litr2c[layer];
		cs->litr3c[layer] -= cf->litr3c_to_litr2c[layer];
		/* Fluxes out of lignin litter pool */
		cs->HRlitr4_snk += cf->litr4_hr[layer];
		cs->litr4c[layer] -= cf->litr4_hr[layer];
		cs->soil3c[layer] += cf->litr4c_to_soil3c[layer];
		cs->litr4c[layer] -= cf->litr4c_to_soil3c[layer];
		/* Fluxes out of fast soil pool */
		cs->HRsoil1_snk += cf->soil1_hr[layer];
		cs->soil1c[layer] -= cf->soil1_hr[layer];
		cs->soil2c[layer] += cf->soil1c_to_soil2c[layer];
		cs->soil1c[layer] -= cf->soil1c_to_soil2c[layer];
		/* Fluxes out of medium soil pool */
		cs->HRsoil2_snk += cf->soil2_hr[layer];
		cs->soil2c[layer] -= cf->soil2_hr[layer];
		cs->soil3c[layer] += cf->soil2c_to_soil3c[layer];
		cs->soil2c[layer] -= cf->soil2c_to_soil3c[layer];
		/* Fluxes out of slow soil pool */
		cs->HRsoil3_snk += cf->soil3_hr[layer];
		cs->soil3c[layer] -= cf->soil3_hr[layer];
		cs->soil4c[layer] += cf->soil3c_to_soil4c[layer];
		cs->soil3c[layer] -= cf->soil3c_to_soil4c[layer];
		/* Fluxes out of recalcitrant SOM pool */
		cs->HRsoil4_snk += cf->soil4_hr[layer];
		cs->soil4c[layer] -= cf->soil4_hr[layer];

		/* Fluxes out of coarse woody debris into litter pools */
		ns->litr2n[layer] += nf->cwdn_to_litr2n[layer];
		ns->litr3n[layer] += nf->cwdn_to_litr3n[layer];
		ns->litr4n[layer] += nf->cwdn_to_litr4n[layer];
		ns->cwdn[layer] -= (nf->cwdn_to_litr2n[layer] + nf->cwdn_to_litr3n[layer] + nf->cwdn_to_litr4n[layer]);
		/* Fluxes out of labile litter pool */
		ns->soil1n[layer] += nf->litr1n_to_soil1n[layer];
		ns->litr1n[layer] -= nf->litr1n_to_soil1n[layer];
		ns->litr1n[layer] -= nf->litr1n_to_release[layer];
		/* Fluxes out of cellulose litter pool */
		ns->soil2n[layer] += nf->litr2n_to_soil2n[layer];
		ns->litr2n[layer] -= nf->litr2n_to_soil2n[layer];
		ns->litr2n[layer] -= nf->litr2n_to_release[layer];
		/* Fluxes from shielded to unshielded cellulose pools */
		ns->litr2n[layer] += nf->litr3n_to_litr2n[layer];
		ns->litr3n[layer] -= nf->litr3n_to_litr2n[layer];
		/* Fluxes out of lignin litter pool */
		ns->soil3n[layer] += nf->litr4n_to_soil3n[layer];
		ns->litr4n[layer] -= nf->litr4n_to_soil3n[layer];
		ns->litr4n[layer] -= nf->litr4n_to_release[layer];
		/* Fluxes out of fast soil pool */
		ns->soil2n[layer] += nf->soil1n_to_soil2n[layer];
		ns->soil1n[layer] -= nf->soil1n_to_soil2n[layer];
		/* Fluxes out of medium soil pool */
		ns->soil3n[layer] += nf->soil2n_to_soil3n[layer];
		ns->soil2n[layer] -= nf->soil2n_to_soil3n[layer];
		/* Fluxes out of slow soil pool */
		ns->soil4n[layer] += nf->soil3n_to_soil4n[layer];
		ns->soil3n[layer] -= nf->soil3n_to_soil4n[layer];

		/*-------------------------------------------------*/
		/* 8.1 control */
		if ((ns->soil1n[layer] == 0 && fabs(cs->soil1c[layer]) > CRIT_PREC_lenient) || (ns->soil2n[layer] == 0 && fabs(cs->soil2c[layer]) > CRIT_PREC_lenient) ||
			(ns->soil3n[layer] == 0 && fabs(cs->soil3c[layer]) > CRIT_PREC_lenient) || (ns->soil4n[layer] == 0 && fabs(cs->soil4c[layer]) > CRIT_PREC_lenient) ||
			(ns->litr1n[layer] == 0 && fabs(cs->litr1c[layer]) > CRIT_PREC_lenient) || (ns->litr2n[layer] == 0 && fabs(cs->litr2c[layer]) > CRIT_PREC_lenient) ||
			(ns->litr3n[layer] == 0 && fabs(cs->litr3c[layer]) > CRIT_PREC_lenient) || (ns->litr4n[layer] == 0 && fabs(cs->litr4c[layer]) > CRIT_PREC_lenient))
		{
			if (!errorCode)
			{
				printf("\n");
				printf("ERROR: CN ratio problem in state_update.c\n");
				errorCode = 1;
			}
		}
		/*-------------------------------------------------*/
		/* 8.2 soil pool changes for GW-calculations */
	
	
		soilInfo->dismatTOTALdecomp[2][layer] = nf->litr1n_to_soil1n[layer] - nf->soil1n_to_soil2n[layer];
		soilInfo->dismatTOTALdecomp[3][layer] = nf->litr2n_to_soil2n[layer] + nf->soil1n_to_soil2n[layer] - nf->soil2n_to_soil3n[layer];
		soilInfo->dismatTOTALdecomp[4][layer] = nf->litr4n_to_soil3n[layer] + nf->soil2n_to_soil3n[layer] - nf->soil3n_to_soil4n[layer];
		soilInfo->dismatTOTALdecomp[5][layer] = nf->soil3n_to_soil4n[layer];
		soilInfo->dismatTOTALdecomp[6][layer] = cf->litr1c_to_soil1c[layer] - cf->soil1c_to_soil2c[layer] - cf->soil1_hr[layer];
		soilInfo->dismatTOTALdecomp[7][layer] = cf->litr2c_to_soil2c[layer] + cf->soil1c_to_soil2c[layer] - cf->soil2c_to_soil3c[layer] - cf->soil2_hr[layer];
		soilInfo->dismatTOTALdecomp[8][layer] = cf->litr4c_to_soil3c[layer] + cf->soil2c_to_soil3c[layer] - cf->soil3c_to_soil4c[layer] - cf->soil3_hr[layer];
		soilInfo->dismatTOTALdecomp[9][layer] = cf->soil3c_to_soil4c[layer] - cf->soil4_hr[layer];


		if (layer < GWlayer)
		{
			soilInfo->dismatUNSATdecomp[2][layer] = soilInfo->dismatTOTALdecomp[2][layer];
			soilInfo->dismatUNSATdecomp[3][layer] = soilInfo->dismatTOTALdecomp[3][layer];
			soilInfo->dismatUNSATdecomp[4][layer] = soilInfo->dismatTOTALdecomp[4][layer];
			soilInfo->dismatUNSATdecomp[5][layer] = soilInfo->dismatTOTALdecomp[5][layer];
			soilInfo->dismatUNSATdecomp[6][layer] = soilInfo->dismatTOTALdecomp[6][layer];
			soilInfo->dismatUNSATdecomp[7][layer] = soilInfo->dismatTOTALdecomp[7][layer];
			soilInfo->dismatUNSATdecomp[8][layer] = soilInfo->dismatTOTALdecomp[8][layer];
			soilInfo->dismatUNSATdecomp[9][layer] = soilInfo->dismatTOTALdecomp[9][layer];

		}

		/* transfer value: content_array -> NH4, NO3, DOC, DON */
		if (!errorCode && check_soilcontent(layer, 0, sprop, cs, ns, soilInfo))
		{
			printf("ERROR in check_soilcontent.c for state_update.c\n");
			errorCode = 1;
		}
		
	}



	/*-------------------------------------------------*/
	/* 8.4. Groundwater transport */

	if (sprop->GWlayer != DATA_GAP)
	{ 

		/* if capillary zone exists in unsaturated zone (not in GWlayer) */
		if (sprop->dz_CAPILcf+sprop->dz_NORMcf)
		{
			for (dm = 0; dm < N_DISSOLVMATER; dm++)
			{
				/* calculation of NORM and CAPIL ratio */
				if (soilInfo->content_NORMcf[dm] + soilInfo->content_CAPILcf[dm])
				{
					ratioNORM = soilInfo->content_NORMcf[dm] / (soilInfo->content_NORMcf[dm] + soilInfo->content_CAPILcf[dm]);
					ratioCAPIL = soilInfo->content_CAPILcf[dm] / (soilInfo->content_NORMcf[dm] + soilInfo->content_CAPILcf[dm]);
				}
				else
				{
					ratioNORM = sprop->dz_NORMcf / (sprop->dz_NORMcf + sprop->dz_CAPILcf);
					ratioCAPIL = sprop->dz_CAPILcf / (sprop->dz_NORMcf + sprop->dz_CAPILcf);
				}
				if (fabs(1 - ratioNORM - ratioCAPIL) > CRIT_PREC)
				{
					printf("\n");
					printf("ERROR in ratio calculation in capillary_diffusion.c\n");
					errorCode = 1;
				}

				diffNORM = soilInfo->dismatTOTALdecomp[dm][CFlayer] * ratioNORM;
				diffCAPIL = soilInfo->dismatTOTALdecomp[dm][CFlayer] * ratioCAPIL;

				/* NORM zone: negative storage value is not possible  */
				if (soilInfo->content_NORMcf[dm] + diffNORM < 0) diffNORM = -1 * soilInfo->content_NORMcf[dm];

				soilInfo->content_NORMcf[dm] += diffNORM;

				diffCAPIL = soilInfo->dismatTOTALdecomp[dm][CFlayer] - diffNORM;

				/* CAPIL zone: negative storage value is not possible */
				if (soilInfo->content_CAPILcf[dm] + diffCAPIL < 0)
				{
					if (fabs(soilInfo->content_CAPILcf[dm] + diffCAPIL) > CRIT_PREC)
					{
						printf("ERROR in content_CAPILcf calculation in state_update.c\n");
						errorCode = 1;
					}
					else
						diffCAPIL = soilInfo->content_CAPILcf[dm];
				}
				soilInfo->content_CAPILcf[dm] += diffCAPIL;

				if (fabs(soilInfo->content_soil[dm][CFlayer] - (soilInfo->content_NORMcf[dm] + soilInfo->content_CAPILcf[dm])) > CRIT_PREC_lenient*10)
				{
					printf("ERROR in content_CAPILcf calculation in state_update.c\n");
			        errorCode = 1;
				}
		
			}

		}
		for (layer = GWlayer; layer < N_SOILLAYERS; layer++)
		{
			/* in GWlayer: soilInfo->dismatTOTALdecomp is distributed proportionally - soilInfo->dismatTOTALdecomp in saturated zone is accounted for  GWdecomp: -: sink, +: source*/

			if (layer == GWlayer)
			{
	
				for (dm = 0; dm < N_DISSOLVMATER; dm++)
				{
					if (soilInfo->content_NORMgw[dm] + soilInfo->content_CAPILgw[dm] + soilInfo->content_SATgw[dm])
					{
						ratioNORM = soilInfo->content_NORMgw[dm] / (soilInfo->content_NORMgw[dm] + soilInfo->content_CAPILgw[dm] + soilInfo->content_SATgw[dm]);
						ratioCAPIL = soilInfo->content_CAPILgw[dm] / (soilInfo->content_NORMgw[dm] + soilInfo->content_CAPILgw[dm] + soilInfo->content_SATgw[dm]);
						ratioSAT = soilInfo->content_SATgw[dm] / (soilInfo->content_NORMgw[dm] + soilInfo->content_CAPILgw[dm] + soilInfo->content_SATgw[dm]);
					}
					else
					{
						ratioNORM = sprop->dz_NORMgw / sitec->soillayer_thickness[GWlayer];
						ratioCAPIL = sprop->dz_CAPILgw / sitec->soillayer_thickness[GWlayer];
						ratioSAT = sprop->dz_SATgw / sitec->soillayer_thickness[GWlayer];
					}
					if (fabs(1 - ratioNORM - ratioCAPIL - ratioSAT) > CRIT_PREC)
					{
						printf("\n");
						printf("ERROR in ratio calculation in capillary_diffusion.c\n");
						errorCode = 1;
					}

					if (soilInfo->dismatTOTALdecomp[dm][layer])
					{
						diffNORM = soilInfo->dismatTOTALdecomp[dm][layer] * ratioNORM;
						diffCAPIL = soilInfo->dismatTOTALdecomp[dm][layer] * ratioCAPIL;
						diffSAT = soilInfo->dismatTOTALdecomp[dm][layer] * ratioSAT;


						/* NORM zone (negative content is not possible - covered by groundwater) */
						soilInfo->contentBOUND_NORMgw[dm] += diffNORM * (1 - soilInfo->dissolv_prop[dm]);
						soilInfo->contentDISSOLV_NORMgw[dm] += diffNORM * soilInfo->dissolv_prop[dm];
						if (soilInfo->contentBOUND_NORMgw[dm] < 0)
						{
							diffNORM -= soilInfo->contentBOUND_NORMgw[dm];
							soilInfo->contentBOUND_NORMgw[dm] = 0;
						}
						if (soilInfo->contentDISSOLV_NORMgw[dm] < 0)
						{
							diffNORM -= soilInfo->contentDISSOLV_NORMgw[dm];
							soilInfo->contentDISSOLV_NORMgw[dm] = 0;
						}

						/* CAPIL zone (negative content is not possible - covered by groundwater) */
						soilInfo->contentBOUND_CAPILgw[dm] += diffCAPIL * (1 - soilInfo->dissolv_prop[dm]);
						soilInfo->contentDISSOLV_CAPILgw[dm] += diffCAPIL * soilInfo->dissolv_prop[dm];
						if (soilInfo->contentBOUND_CAPILgw[dm] < 0)
						{
							diffCAPIL -= soilInfo->contentBOUND_CAPILgw[dm];
							soilInfo->contentBOUND_CAPILgw[dm] = 0;
						}
						if (soilInfo->contentDISSOLV_CAPILgw[dm] < 0)
						{
							diffCAPIL -= soilInfo->contentDISSOLV_CAPILgw[dm];
							soilInfo->contentDISSOLV_CAPILgw[dm] = 0;
						}

						/* SAT zone: diffSAT is the part remained after diffNORM and diffCAPIL, but: dissolved part of SATgw is unchanged - GWdecomp */
						diffSAT = soilInfo->dismatTOTALdecomp[dm][layer] - diffNORM - diffCAPIL;
						
						soilInfo->contentBOUND_SATgw[dm] += diffSAT * (1 - soilInfo->dissolv_prop[dm]);
						soilInfo->dismatGWdecomp[dm][GWlayer] = -1 * diffSAT * soilInfo->dissolv_prop[dm];

						/* negative content is not possible - covered by groundwater */
						if (soilInfo->contentBOUND_SATgw[dm] < 0)
						{
							soilInfo->dismatGWdecomp[dm][layer] -= soilInfo->contentBOUND_SATgw[dm];
							soilInfo->contentBOUND_SATgw[dm] = 0;
						}

						/* change of UNSAT zone*/
						soilInfo->dismatUNSATdecomp[dm][layer] = diffNORM + diffCAPIL;

						
					}
				}

			}
			else
			{
				/* below GWlayer - dissolv part of dismatTOTALdecomp is covered by GW*/
				for (dm = 0; dm < N_DISSOLVMATER; dm++)
				{
					soilInfo->contentBOUND_soil[dm][layer] += soilInfo->dismatTOTALdecomp[dm][layer] * (1 - soilInfo->dissolv_prop[dm]);
					soilInfo->dismatGWdecomp[dm][layer] = -1 * soilInfo->dismatTOTALdecomp[dm][layer] * soilInfo->dissolv_prop[dm];

					/* negative content is not possible - covered by groundwater */
					if (soilInfo->contentBOUND_soil[dm][layer] < 0)
					{
						soilInfo->dismatGWdecomp[dm][layer] -= soilInfo->contentBOUND_soil[dm][layer];
						soilInfo->contentBOUND_soil[dm][layer] = 0;
					}
				}

			}
		}

	}


	/* *****************************************************************************************************/
	/* src/snk variables*/
	
	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		for (dm = 0; dm < N_DISSOLVN; dm++) soilInfo->dismatGWdecompN_total += soilInfo->dismatGWdecomp[dm][layer];
		for (dm = N_DISSOLVN; dm < N_DISSOLVMATER; dm++) soilInfo->dismatGWdecompC_total += soilInfo->dismatGWdecomp[dm][layer];
	}


	if (soilInfo->dismatGWdecompN_total > 0)
		ns->GWsrc_N += soilInfo->dismatGWdecompN_total;
	else
		ns->GWsnk_N += -1 * soilInfo->dismatGWdecompN_total;


	if (soilInfo->dismatGWdecompC_total > 0)
		cs->GWsrc_C += soilInfo->dismatGWdecompC_total;
	else
		cs->GWsnk_C += -1 * soilInfo->dismatGWdecompC_total;


	/* *****************************************************************************************************/
	/* update pools */


	/*  firsttime_flag=1 (after calculation note initial values, int partlyORtotal_flag=1 (TOTAL (BOUND+DISSOLV) is affected  */
	if (!errorCode && calc_DISSOLVandBOUND(1, 1, sprop, soilInfo))
	{
		printf("ERROR in calc_DISSOLVandBOUND.c for state_update.c\n");
		errorCode = 1;
	}

	if (!errorCode && check_soilcontent(-1, 1, sprop, cs, ns, soilInfo))
	{
		printf("ERROR in check_soilcontent.c for state_update.c\n");
		errorCode = 1;
	}



	/***************************************************************************************************************************************************/
	/* 9. Daily allocation fluxes */
	/* daily leaf allocation fluxes */
	
	if (!errorCode && epc->leaf_cn && CNratio_control(cs, epc->leaf_cn, cs->leafc, ns->leafn, cf->cpool_to_leafc, nf->npool_to_leafn, 0)) 
	{
		cs->cpool          -= cf->cpool_to_leafc;
		cs->leafc          += cf->cpool_to_leafc;
		
		ns->npool          -= nf->npool_to_leafn;
		ns->leafn           = cs->leafc / epc->leaf_cn;
	}
	else if (epc->leaf_cn)		
	{
		if (!errorCode) printf("ERROR in cpool_to_leafc CN calculation in state_update.c\n");
		errorCode=1;
	}

	if (!errorCode && epc->leaf_cn && CNratio_control(cs, epc->leaf_cn, cs->leafc_storage, ns->leafn_storage, cf->cpool_to_leafc_storage, nf->npool_to_leafn_storage, 0)) 
	{
		cs->cpool          -= cf->cpool_to_leafc_storage;
		cs->leafc_storage  += cf->cpool_to_leafc_storage;
		
		ns->npool          -= nf->npool_to_leafn_storage;
		ns->leafn_storage   = cs->leafc_storage / epc->leaf_cn;
	}
	else if (epc->leaf_cn)		
	{
		if (!errorCode) printf("ERROR in cpool_to_leafc_storage CN calculation in state_update.c\n");
		errorCode=1;
	}

	/* Daily fine root allocation fluxes */
	if (!errorCode && epc->froot_cn && CNratio_control(cs, epc->froot_cn, cs->frootc, ns->frootn, cf->cpool_to_frootc, nf->npool_to_frootn, 0)) 
	{
		cs->cpool          -= cf->cpool_to_frootc;
		cs->frootc         += cf->cpool_to_frootc;
		
		ns->npool          -= nf->npool_to_frootn;
		ns->frootn          = cs->frootc / epc->froot_cn;
	}
	else if (epc->froot_cn)		
	{
		if (!errorCode) printf("ERROR in cpool_to_frootc CN calculation in state_update.c\n");
		errorCode=1;
	}

	if (!errorCode && epc->froot_cn && CNratio_control(cs, epc->froot_cn, cs->frootc_storage, ns->frootn_storage, cf->cpool_to_frootc_storage, nf->npool_to_frootn_storage, 0)) 
	{
		cs->cpool          -= cf->cpool_to_frootc_storage;
		cs->frootc_storage += cf->cpool_to_frootc_storage;
		
		ns->npool          -= nf->npool_to_frootn_storage;
		ns->frootn_storage  = cs->frootc_storage / epc->froot_cn;
	}
	else if (epc->froot_cn)		
	{
		if (!errorCode) printf("ERROR in cpool_to_frootc_storage CN calculation in state_update.c\n");
		errorCode=1;
	}

	/* Daily yield allocation fluxes + EXTRA: effect of heat stress during flowering in yield filling phenophase */
	if (!errorCode && epc->yield_cn && CNratio_control(cs, epc->yield_cn, cs->yieldc, ns->yieldn, cf->cpool_to_yield, nf->npool_to_yieldn, 0)) 
	{
		cs->cpool          -= cf->cpool_to_yield;
		cs->yieldc         += cf->cpool_to_yield;
		
		ns->npool          -= nf->npool_to_yieldn;
		ns->yieldn          = cs->yieldc / epc->yield_cn;
			
		cs->yieldc         -= cf->yieldc_to_flowHS;
		cs->STDBc_yield    += cf->yieldc_to_flowHS;

		ns->yieldn         -= nf->yieldn_to_flowHS;
		ns->STDBn_yield    += nf->yieldn_to_flowHS;

		/* control */
		if ((cf->yieldc_to_flowHS > 0 && nf->yieldn_to_flowHS > 0 && epv->n_actphen != epc->n_flowHS_phenophase) || 
			(cs->yieldc < 0 && fabs(cs->yieldc) > CRIT_PREC)  || (ns->yieldn < 0 && fabs(ns->yieldn) > CRIT_PREC)  )
		{
			if (!errorCode) printf("ERROR in flowering heat stress calculation in state_update.c\n");
			errorCode=1;
		}
	}
	else if (epc->yield_cn)		
	{
		if (!errorCode) printf("ERROR in cpool_to_yield CN calculation in state_update.c\n");
		errorCode=1;
	}

	if (!errorCode && epc->yield_cn && CNratio_control(cs, epc->yield_cn, cs->yieldc_storage, ns->yieldn_storage, cf->cpool_to_yieldc_storage, nf->npool_to_yieldn_storage, 0)) 
	{
		cs->cpool          -= cf->cpool_to_yieldc_storage;
		cs->yieldc_storage += cf->cpool_to_yieldc_storage;
		
		ns->npool          -= nf->npool_to_yieldn_storage;
		ns->yieldn_storage  = cs->yieldc_storage / epc->yield_cn;
	}
	else if (epc->yield_cn)		
	{
		if (!errorCode) printf("ERROR in cpool_to_yieldc_storage CN calculation in state_update.c\n");
		errorCode=1;
	}
	
	/* Daily sofstem allocation fluxes */
	if (!errorCode && epc->softstem_cn && CNratio_control(cs, epc->softstem_cn, cs->softstemc, ns->softstemn, cf->cpool_to_softstemc, nf->npool_to_softstemn, 0)) 
	{
		cs->cpool          -= cf->cpool_to_softstemc;
		cs->softstemc       += cf->cpool_to_softstemc;
		
		ns->npool          -= nf->npool_to_softstemn;
		ns->softstemn       = cs->softstemc / epc->softstem_cn;
	}
	else if (epc->softstem_cn)		
	{
		if (!errorCode) printf("ERROR in cpool_to_softstemc CN calculation in state_update.c\n");
		errorCode=1;
	}

	if (!errorCode && epc->softstem_cn && CNratio_control(cs, epc->softstem_cn, cs->softstemc_storage, ns->softstemn_storage, cf->cpool_to_softstemc_storage, nf->npool_to_softstemn_storage, 0)) 
	{
		cs->cpool              -= cf->cpool_to_softstemc_storage;
		cs->softstemc_storage  += cf->cpool_to_softstemc_storage;
		
		ns->npool              -= nf->npool_to_softstemn_storage;
		ns->softstemn_storage   = cs->softstemc_storage / epc->softstem_cn;
	}
	else if (epc->softstem_cn)		
	{
		if (!errorCode) printf("ERROR in cpool_to_softstemc_storage CN calculation in state_update.c\n");
		errorCode=1;
	}			

	/* Daily live stem wood allocation fluxes */
	if (!errorCode && epc->livewood_cn && CNratio_control(cs, epc->livewood_cn, cs->livestemc, ns->livestemn, cf->cpool_to_livestemc, nf->npool_to_livestemn, 0)) 
	{
		cs->cpool      -= cf->cpool_to_livestemc;
		cs->livestemc  += cf->cpool_to_livestemc;
		
		ns->npool      -= nf->npool_to_livestemn;
		ns->livestemn   = cs->livestemc / epc->livewood_cn;
	}
	else if (epc->livewood_cn)		
	{
		if (!errorCode) printf("ERROR in cpool_to_livestemc CN calculation in state_update.c\n");
		errorCode=1;
	}

	if (!errorCode && epc->livewood_cn && CNratio_control(cs, epc->livewood_cn, cs->livestemc_storage, ns->livestemn_storage, cf->cpool_to_livestemc_storage, nf->npool_to_livestemn_storage, 0)) 
	{
		cs->cpool              -= cf->cpool_to_livestemc_storage;
		cs->livestemc_storage  += cf->cpool_to_livestemc_storage;
		
		ns->npool              -= nf->npool_to_livestemn_storage;
		ns->livestemn_storage   = cs->livestemc_storage / epc->livewood_cn;
	}
	else if (epc->livewood_cn)		
	{
		if (!errorCode) printf("ERROR in cpool_to_livestemc_storage CN calculation in state_update.c\n");
		errorCode=1;
	}
	
	/* Daily dead stem wood allocation fluxes */
	if (!errorCode && epc->deadwood_cn && CNratio_control(cs, epc->deadwood_cn, cs->deadstemc, ns->deadstemn, cf->cpool_to_deadstemc, nf->npool_to_deadstemn, 0)) 
	{
		cs->cpool          -= cf->cpool_to_deadstemc;
		cs->deadstemc          += cf->cpool_to_deadstemc;
		
		ns->npool          -= nf->npool_to_deadstemn;
		ns->deadstemn           = cs->deadstemc / epc->deadwood_cn;
	}
	else if (epc->deadwood_cn)		
	{
		if (!errorCode) printf("ERROR in cpool_to_deadstemc CN calculation in state_update.c\n");
		errorCode=1;
	}

	if (!errorCode && epc->deadwood_cn && CNratio_control(cs, epc->deadwood_cn, cs->deadstemc_storage, ns->deadstemn_storage, cf->cpool_to_deadstemc_storage, nf->npool_to_deadstemn_storage, 0)) 
	{
		cs->cpool          -= cf->cpool_to_deadstemc_storage;
		cs->deadstemc_storage  += cf->cpool_to_deadstemc_storage;
		
		ns->npool          -= nf->npool_to_deadstemn_storage;
		ns->deadstemn_storage   = cs->deadstemc_storage / epc->deadwood_cn;
	}
	else if (epc->deadwood_cn)		
	{
		if (!errorCode) printf("ERROR in cpool_to_deadstemc_storage CN calculation in state_update.c\n");
		errorCode=1;
	}
	
	/* Daily live coarse root wood allocation fluxes */
	if (!errorCode && epc->livewood_cn && CNratio_control(cs, epc->livewood_cn, cs->livecrootc, ns->livecrootn, cf->cpool_to_livecrootc, nf->npool_to_livecrootn, 0)) 
	{
		cs->cpool          -= cf->cpool_to_livecrootc;
		cs->livecrootc          += cf->cpool_to_livecrootc;
		
		ns->npool          -= nf->npool_to_livecrootn;
		ns->livecrootn           = cs->livecrootc / epc->livewood_cn;
	}
	else if (epc->livewood_cn)		
	{
		if (!errorCode) printf("ERROR in cpool_to_livecrootc CN calculation in state_update.c\n");
		errorCode=1;
	}

	if (!errorCode && epc->livewood_cn && CNratio_control(cs, epc->livewood_cn, cs->livecrootc_storage, ns->livecrootn_storage, cf->cpool_to_livecrootc_storage, nf->npool_to_livecrootn_storage, 0)) 
	{
		cs->cpool          -= cf->cpool_to_livecrootc_storage;
		cs->livecrootc_storage  += cf->cpool_to_livecrootc_storage;
		
		ns->npool          -= nf->npool_to_livecrootn_storage;
		ns->livecrootn_storage   = cs->livecrootc_storage / epc->livewood_cn;
	}
	else if (epc->livewood_cn)		
	{
		if (!errorCode) printf("ERROR in cpool_to_livecrootc_storage CN calculation in state_update.c\n");
		errorCode=1;
	}

	/* Daily dead coarse root wood allocation fluxes */
	if (!errorCode && epc->deadwood_cn && CNratio_control(cs, epc->deadwood_cn, cs->deadcrootc, ns->deadcrootn, cf->cpool_to_deadcrootc, nf->npool_to_deadcrootn, 0)) 
	{
		cs->cpool          -= cf->cpool_to_deadcrootc;
		cs->deadcrootc          += cf->cpool_to_deadcrootc;
		
		ns->npool          -= nf->npool_to_deadcrootn;
		ns->deadcrootn           = cs->deadcrootc / epc->deadwood_cn;
	}
	else if (epc->deadwood_cn)		
	{
		if (!errorCode) printf("ERROR in cpool_to_deadcrootc CN calculation in state_update.c\n");
		errorCode=1;
	}

	if (!errorCode && epc->deadwood_cn && CNratio_control(cs, epc->deadwood_cn, cs->deadcrootc_storage, ns->deadcrootn_storage, cf->cpool_to_deadcrootc_storage, nf->npool_to_deadcrootn_storage, 0)) 
	{
		cs->cpool          -= cf->cpool_to_deadcrootc_storage;
		cs->deadcrootc_storage  += cf->cpool_to_deadcrootc_storage;
		
		ns->npool          -= nf->npool_to_deadcrootn_storage;
		ns->deadcrootn_storage   = cs->deadcrootc_storage / epc->deadwood_cn;
	}
	else if (epc->deadwood_cn)		
	{
		if (!errorCode) printf("ERROR in cpool_to_deadcrootc_storage CN calculation in state_update.c\n");
		errorCode=1;
	}
	

	/* Daily allocation for transfer growth respiration */
	cs->gresp_storage  += cf->cpool_to_gresp_storage;
	cs->cpool          -= cf->cpool_to_gresp_storage;
	 
	/* Calculate sum of leafC in a given phenophase */
	pp = (int) epv->n_actphen - 1;
	cs->leafcSUM_phenphase[pp] += cf->cpool_to_leafc + cf->leafc_transfer_to_leafc;

	/***************************************************************************************************************************************************/
	/* 10. Daily growth respiration fluxes */
	/* Leaf growth respiration */
	cs->GRleaf_snk     += cf->cpool_leaf_GR;
	cs->cpool           -= cf->cpool_leaf_GR;
	cs->GRleaf_snk     += cf->cpool_leaf_storage_GR;
	cs->cpool           -= cf->cpool_leaf_storage_GR;
	cs->GRleaf_snk     += cf->transfer_leaf_GR;
	cs->gresp_transfer  -= cf->transfer_leaf_GR;
	/* Fine root growth respiration */
	cs->GRfroot_snk    += cf->cpool_froot_GR;
	cs->cpool           -= cf->cpool_froot_GR;
	cs->GRfroot_snk    += cf->cpool_froot_storage_GR;
	cs->cpool           -= cf->cpool_froot_storage_GR;
	cs->GRfroot_snk    += cf->transfer_froot_GR;
	cs->gresp_transfer  -= cf->transfer_froot_GR;
	/* yield growth respiration */
	cs->GRyield_snk     += cf->cpool_yield_GR;
	cs->cpool            -= cf->cpool_yield_GR;
	cs->GRyield_snk     += cf->cpool_yieldc_storage_GR;
	cs->cpool            -= cf->cpool_yieldc_storage_GR;
	cs->GRyield_snk     += cf->transfer_yield_GR;
	cs->gresp_transfer   -= cf->transfer_yield_GR;
	/* yield growth respiration. */
	cs->GRsoftstem_snk  += cf->cpool_softstem_GR;
	cs->cpool            -= cf->cpool_softstem_GR;
	cs->GRsoftstem_snk  += cf->cpool_softstem_storage_GR;
	cs->cpool            -= cf->cpool_softstem_storage_GR;
	cs->GRsoftstem_snk  += cf->transfer_softstem_GR;
	cs->gresp_transfer   -= cf->transfer_softstem_GR;
	/* Live stem growth respiration */ 
	cs->GRlivestem_snk  += cf->cpool_livestem_GR;
	cs->cpool            -= cf->cpool_livestem_GR;
	cs->GRlivestem_snk  += cf->cpool_livestem_storage_GR;
	cs->cpool            -= cf->cpool_livestem_storage_GR;
	cs->GRlivestem_snk  += cf->transfer_livestem_GR;
	cs->gresp_transfer   -= cf->transfer_livestem_GR;
	/* Dead stem growth respiration */ 
	cs->GRdeadstem_snk  += cf->cpool_deadstem_GR;
	cs->cpool            -= cf->cpool_deadstem_GR;
	cs->GRdeadstem_snk  += cf->cpool_deadstem_storage_GR;
	cs->cpool            -= cf->cpool_deadstem_storage_GR;
	cs->GRdeadstem_snk  += cf->transfer_deadstem_GR;
	cs->gresp_transfer   -= cf->transfer_deadstem_GR;
	/* Live coarse root growth respiration */ 
	cs->GRlivecroot_snk += cf->cpool_livecroot_GR;
	cs->cpool            -= cf->cpool_livecroot_GR;
	cs->GRlivecroot_snk += cf->cpool_livecroot_storage_GR;
	cs->cpool            -= cf->cpool_livecroot_storage_GR;
	cs->GRlivecroot_snk += cf->transfer_livecroot_GR;
	cs->gresp_transfer   -= cf->transfer_livecroot_GR;
	/* Dead coarse root growth respiration */ 
	cs->GRdeadcroot_snk += cf->cpool_deadcroot_GR;
	cs->cpool            -= cf->cpool_deadcroot_GR;
	cs->GRdeadcroot_snk += cf->cpool_deadcroot_storage_GR;
	cs->cpool            -= cf->cpool_deadcroot_storage_GR;
	cs->GRdeadcroot_snk += cf->transfer_deadcroot_GR;
	cs->gresp_transfer   -= cf->transfer_deadcroot_GR;

	
	/***************************************************************************************************************************************************/
	/* 11. Maintanance respiration fluxes
	       covering of maintananance respiration fluxes from NSC (non-structural carbohydrate), namely storage and transfer pools*/
	if (!errorCode && MRdeficit_calculation(epc, ctrl, cf, nf, cs, ns))
	{
		errorCode=1;
		if (!errorCode) printf("ERROR in MRdeficit_calculation.c from bgc.c\n");
	}
	
	
	/***************************************************************************************************************************************************/
	
	/* 12. Annual allocation fluxes, one day per year */
	if (alloc)
	{
		/* Move storage material into transfer compartments on the annual
		allocation day. This is a special case, where a flux is assessed in
		the state_variable update routine.  This is required to have the
		allocation of excess C and N show up as new growth in the next growing
		season, instead of two growing seasons from now. */
		cf->leafc_storage_to_leafc_transfer				= cs->leafc_storage;
		cf->frootc_storage_to_frootc_transfer			= cs->frootc_storage;
		cf->yieldc_storage_to_yieldc_transfer			= cs->yieldc_storage;
		cf->softstemc_storage_to_softstemc_transfer		= cs->softstemc_storage;
		cf->gresp_storage_to_gresp_transfer				= cs->gresp_storage;
		cf->livestemc_storage_to_livestemc_transfer		= cs->livestemc_storage;
		cf->deadstemc_storage_to_deadstemc_transfer		= cs->deadstemc_storage;
		cf->livecrootc_storage_to_livecrootc_transfer	= cs->livecrootc_storage;
		cf->deadcrootc_storage_to_deadcrootc_transfer	= cs->deadcrootc_storage;

		nf->leafn_storage_to_leafn_transfer				= ns->leafn_storage;
		nf->frootn_storage_to_frootn_transfer			= ns->frootn_storage;
		nf->yieldn_storage_to_yieldn_transfer			= ns->yieldn_storage;
		nf->softstemn_storage_to_softstemn_transfer		= ns->softstemn_storage;
		nf->livestemn_storage_to_livestemn_transfer		= ns->livestemn_storage;
		nf->deadstemn_storage_to_deadstemn_transfer		= ns->deadstemn_storage;
		nf->livecrootn_storage_to_livecrootn_transfer	= ns->livecrootn_storage;
		nf->deadcrootn_storage_to_deadcrootn_transfer	= ns->deadcrootn_storage;
		

		/* update states variables */
		cs->leafc_transfer     += cf->leafc_storage_to_leafc_transfer;
		cs->leafc_storage      -= cf->leafc_storage_to_leafc_transfer;
		cs->frootc_transfer    += cf->frootc_storage_to_frootc_transfer;
		cs->frootc_storage     -= cf->frootc_storage_to_frootc_transfer;
		cs->yieldc_transfer    += cf->yieldc_storage_to_yieldc_transfer;
		cs->yieldc_storage     -= cf->yieldc_storage_to_yieldc_transfer;
		cs->softstemc_transfer += cf->softstemc_storage_to_softstemc_transfer;
		cs->softstemc_storage  -= cf->softstemc_storage_to_softstemc_transfer;
		cs->gresp_transfer     += cf->gresp_storage_to_gresp_transfer;
		cs->gresp_storage      -= cf->gresp_storage_to_gresp_transfer;

		ns->leafn_transfer     += nf->leafn_storage_to_leafn_transfer;
		ns->leafn_storage      -= nf->leafn_storage_to_leafn_transfer;
		ns->frootn_transfer    += nf->frootn_storage_to_frootn_transfer;
		ns->frootn_storage     -= nf->frootn_storage_to_frootn_transfer;
		ns->yieldn_transfer    += nf->yieldn_storage_to_yieldn_transfer;
		ns->yieldn_storage     -= nf->yieldn_storage_to_yieldn_transfer;
		ns->softstemn_transfer += nf->softstemn_storage_to_softstemn_transfer;
		ns->softstemn_storage  -= nf->softstemn_storage_to_softstemn_transfer;


		cs->livestemc_transfer  += cf->livestemc_storage_to_livestemc_transfer;
		cs->livestemc_storage   -= cf->livestemc_storage_to_livestemc_transfer;
		cs->deadstemc_transfer  += cf->deadstemc_storage_to_deadstemc_transfer;
		cs->deadstemc_storage   -= cf->deadstemc_storage_to_deadstemc_transfer;
		cs->livecrootc_transfer += cf->livecrootc_storage_to_livecrootc_transfer;
		cs->livecrootc_storage  -= cf->livecrootc_storage_to_livecrootc_transfer;
		cs->deadcrootc_transfer += cf->deadcrootc_storage_to_deadcrootc_transfer;
		cs->deadcrootc_storage  -= cf->deadcrootc_storage_to_deadcrootc_transfer;

		ns->livestemn_transfer  += nf->livestemn_storage_to_livestemn_transfer;
		ns->livestemn_storage   -= nf->livestemn_storage_to_livestemn_transfer;
		ns->deadstemn_transfer  += nf->deadstemn_storage_to_deadstemn_transfer;
		ns->deadstemn_storage   -= nf->deadstemn_storage_to_deadstemn_transfer;
		ns->livecrootn_transfer += nf->livecrootn_storage_to_livecrootn_transfer;
		ns->livecrootn_storage  -= nf->livecrootn_storage_to_livecrootn_transfer;
		ns->deadcrootn_transfer += nf->deadcrootn_storage_to_deadcrootn_transfer;
		ns->deadcrootn_storage  -= nf->deadcrootn_storage_to_deadcrootn_transfer;
		


		/* for deciduous system, force leafc and frootc to exactly 0.0 on the last day */
		if (!evergreen)
		{
			if (-cs->leafc > CRIT_PREC || -cs->frootc > CRIT_PREC || -cs->yieldc > CRIT_PREC || -cs->softstemc > CRIT_PREC)
			{
				printf("\n");
				if (!errorCode) printf("ERROR in state_update.c: negative plant carbon pool\n");
				errorCode=1;
			}
			if (-ns->leafn > CRIT_PREC || -ns->frootn > CRIT_PREC || -ns->yieldn > CRIT_PREC || -ns->softstemn > CRIT_PREC)
			{
				printf("\n");
				if (!errorCode) printf("ERROR in state_update.c: negative plant nitrogen pool\n");
				errorCode=1;
			}
			if (cs->leafc < CRIT_PREC || ns->leafn  < CRIT_PREC) 
			{
				cs->leafc = 0.0;
				ns->leafn = 0.0;
			}
			if (cs->frootc < CRIT_PREC || ns->frootn  < CRIT_PREC) 
			{
				cs->frootc = 0.0;
				ns->frootn = 0.0;
			}
			if (cs->yieldc < CRIT_PREC || ns->yieldn  < CRIT_PREC) 
			{
				cs->yieldc = 0.0;
				ns->yieldn = 0.0;
			}
			if (cs->softstemc < CRIT_PREC || ns->softstemn  < CRIT_PREC) 
			{
				cs->softstemc = 0.0;
				ns->softstemn = 0.0;
			}
		
		}
	} /* end if allocation day */




	return (errorCode);
}			

int CNratio_control(cstate_struct* cs, double CNratio, double cpool, double npool, double cflux, double nflux, double CNratio_flux)
{
	int errorCode = 0;
	double CNdiff = 0;

	if (CNratio_flux == 0) CNratio_flux = CNratio;

	/* control for leaf C:N ratio of pools */
	if ((npool ==0 && cpool > CRIT_PREC ) || (npool > CRIT_PREC  && cpool == 0))
	{
		if (!errorCode) printf("ERROR in CNratio_control for CN_state_update.c\n");
		errorCode = 1;
	}
	
	if(npool > 0 && cpool > 0) CNdiff = cpool/npool - CNratio;
	

	if(fabs(CNdiff) > 0)
	{
		npool = cpool / CNratio;
		if (fabs(CNdiff) > cs->CNratioERR) 
		{
			cs->CNratioERR = fabs(CNdiff);
		}
	
	
	}

	CNdiff=0;

	/* control for leaf C:N ratio of fluxes */
	if ((nflux == 0 && cflux > CRIT_PREC ) || (nflux > CRIT_PREC  && cflux == 0))
	{
		if (!errorCode) printf("ERROR in CNratio_control for CN_state_update.c\n");
		errorCode = 1;
	}

	if(nflux > 0 && cflux > 0) CNdiff = cflux/nflux - CNratio_flux;


	if(fabs(CNdiff) > 0)
	{
		nflux = cflux / CNratio_flux;
		if (fabs(CNdiff) > cs->CNratioERR) cs->CNratioERR = fabs(CNdiff);
		
		CNdiff = cflux/nflux - CNratio;

	}

	return (!errorCode);
}

int MRdeficit_calculation(const epconst_struct* epc, control_struct* ctrl, cflux_struct* cf, nflux_struct* nf, cstate_struct* cs, nstate_struct* ns)
{	
	/* Covering of maintananance respiration fluxes from storage pools */
	
	int errorCode=0;
	double MR, MR_nw, MR_w, NSC, NSCnw, NSCw, SC, SCnw, SCw, NSC_crit, diffTOTAL, diffTOTAL_nw, diffTOTAL_w, diffNSC, diffSC, diff, excess, NSC_avail;
	double reduction, npool_to_MRdef, CpoolNpool_ratio;

	diffTOTAL_nw = diffTOTAL_w = excess = diffSC = diffNSC = diff = reduction = npool_to_MRdef = CpoolNpool_ratio = 0;
	
	/* summarizing maint.resp fluxes and available non-structural carbohydrate fluxes - non-woody and woody */
	MR_nw = cf->leaf_day_MR + cf->leaf_night_MR + cf->froot_MR + cf->yield_MR + cf->softstem_MR;
	MR_w  = cf->livestem_MR + cf->livecroot_MR;
	MR    = MR_nw + MR_w;
	
		
	NSCnw  = (cs->leafc_storage      +  cs->frootc_storage     + cs->yieldc_storage     + cs->softstemc_storage + 
			   cs->leafc_transfer     +  cs->frootc_transfer    + cs->yieldc_transfer    + cs->softstemc_transfer);

	NSCw  = (cs->livestemc_storage  + cs->livecrootc_storage  + cs->deadstemc_storage  + cs->deadcrootc_storage +
		      cs->livestemc_transfer + cs->livecrootc_transfer + cs->deadstemc_transfer + cs->deadcrootc_transfer);

	SCnw = cs->leafc     +  cs->frootc    + cs->yieldc     + cs->softstemc;

	SCw = cs->livestemc  + cs->livecrootc  + cs->deadstemc  + cs->deadcrootc;

	if (fabs(NSCnw) < CRIT_PREC) NSCnw = 0;
	if (fabs(NSCw)  < CRIT_PREC) NSCw = 0;
	if (fabs(SCnw)  < CRIT_PREC) SCnw = 0;
	if (fabs(SCw)   < CRIT_PREC) SCw = 0;

	NSC = NSCnw + NSCw;
	SC  = SCnw  + SCw;

	
	/* calculation of difference between between the demand (MR) and the source (cpool) - non-woody and woody */
	diffTOTAL  = MR - cs->cpool;

	
	/********************************************************************************************************/
	/* I: method 0.: non-separated MR */
	/*******************************************************************************************************/

	if (ctrl->MRdeficit_flag == 0 && MR)
	{

		/* calculation the difference between NSC and diff (based on available amount) */
		if (diffTOTAL > 0)
		{
			/* critical NSC value: NSC pool = fixed ratio of theroretical maximum of NSC value -> but not all is available */
			NSC_crit = (epc->prop_NSCcrit) * (SC * epc->prop_NSCvsSC);
			NSC_avail = NSC - NSC_crit;

			if (NSC_avail < 0) NSC_avail = 0;

			/* calculation of deficit regarding to structural biomass  */
			if (diffTOTAL <= NSC_avail)
			{
				diffNSC = diffTOTAL;
				diffSC = 0;
			}
			else
			{
				diffNSC = NSC_avail;
				diffSC = diffTOTAL - NSC_avail;
			}

			/* calculation of flxues from NSC pools */
			if (NSC)
			{
				cf->leafc_storage_to_MRdef = diffNSC * cs->leafc_storage / NSC;
				cf->frootc_storage_to_MRdef = diffNSC * cs->frootc_storage / NSC;
				cf->yieldc_storage_to_MRdef = diffNSC * cs->yieldc_storage / NSC;
				cf->softstemc_storage_to_MRdef = diffNSC * cs->softstemc_storage / NSC;

				cf->leafc_transfer_to_MRdef = diffNSC * cs->leafc_transfer / NSC;
				cf->frootc_transfer_to_MRdef = diffNSC * cs->frootc_transfer / NSC;
				cf->yieldc_transfer_to_MRdef = diffNSC * cs->yieldc_transfer / NSC;
				cf->softstemc_transfer_to_MRdef = diffNSC * cs->softstemc_transfer / NSC;

				cf->livestemc_storage_to_MRdef = diffNSC * cs->livestemc_storage / NSC;
				cf->livecrootc_storage_to_MRdef = diffNSC * cs->livecrootc_storage / NSC;
				cf->deadstemc_storage_to_MRdef = diffNSC * cs->deadstemc_storage / NSC;
				cf->deadcrootc_storage_to_MRdef = diffNSC * cs->deadcrootc_storage / NSC;

				cf->livestemc_transfer_to_MRdef = diffNSC * cs->livestemc_transfer / NSC;
				cf->livecrootc_transfer_to_MRdef = diffNSC * cs->livecrootc_transfer / NSC;
				cf->deadstemc_transfer_to_MRdef = diffNSC * cs->deadstemc_transfer / NSC;
				cf->deadcrootc_transfer_to_MRdef = diffNSC * cs->deadcrootc_transfer / NSC;
			}
			else
			{
				if (diffNSC)
				{
					printf("ERROR in MRdeficit calculaton in state_update.c\n");
					errorCode = 1;
				}
			}

			/* if NSC is not enough -> transfer from actual pool */
			if (diffSC > 0)
			{
				cf->leafc_to_MRdef = diffSC * ((cf->leaf_day_MR + cf->leaf_night_MR) / MR);
				if (cf->leafc_to_MRdef > cs->leafc)
				{
					cf->leafc_to_MRdef = cs->leafc;
					excess += diffSC * ((cf->leaf_day_MR + cf->leaf_night_MR) / MR) - cf->leafc_to_MRdef;
				}


				cf->frootc_to_MRdef = diffSC * (cf->froot_MR / MR);
				if (cf->frootc_to_MRdef > cs->frootc)
				{
					cf->frootc_to_MRdef = cs->frootc;
					excess += diffSC * (cf->froot_MR / MR) - cf->frootc_to_MRdef;
				}


				cf->yieldc_to_MRdef = diffSC * (cf->yield_MR / MR);
				if (cf->yieldc_to_MRdef > cs->yieldc)
				{
					cf->yieldc_to_MRdef = cs->yieldc;
					excess += diffSC * (cf->yield_MR / MR) - cf->yieldc_to_MRdef;
				}


				cf->softstemc_to_MRdef = diffSC * (cf->softstem_MR / MR);
				if (cf->softstemc_to_MRdef > cs->softstemc)
				{
					cf->softstemc_to_MRdef = cs->softstemc;
					excess += diffSC * (cf->softstem_MR / MR) - cf->softstemc_to_MRdef;
				}

				cf->livestemc_to_MRdef = diffSC * (cf->livestem_MR / MR);
				if (cs->livestemc < cf->livestemc_to_MRdef)
				{
					cf->livestemc_to_MRdef = cs->livestemc;
					excess += diffSC * (cf->livestem_MR / MR) - cf->livestemc_to_MRdef;
				}

				cf->livecrootc_to_MRdef = diffSC * (cf->livecroot_MR / MR);
				if (cs->livecrootc < cf->livecrootc_to_MRdef)
				{
					cf->livecrootc_to_MRdef = cs->livecrootc;
					excess += diffSC * (cf->livecroot_MR / MR) - cf->livecrootc_to_MRdef;
				}


				/* if maintresp of non-woody biomass can not be covered from non-woody biomass -> decreasing of maintanance flux */
				if (excess)
				{
					if (!ctrl->limitMR_flag) ctrl->limitMR_flag = 1;
					cs->MRlimitSUM += excess;
					cf->leaf_day_MR -= excess * (cf->leaf_day_MR / MR);
					cf->leaf_night_MR -= excess * (cf->leaf_night_MR / MR);
					cf->froot_MR -= excess * (cf->froot_MR / MR);
					cf->softstem_MR -= excess * (cf->softstem_MR / MR);
					cf->yield_MR -= excess * (cf->yield_MR / MR);
					cf->livestem_MR -= excess * (cf->livestem_MR / MR);
					cf->livecroot_MR -= excess * (cf->livecroot_MR / MR);
					excess = 0;
				}
			}
		}
	}

	/********************************************************************************************************/
	/* II: method 1.: organ-separated MR */
	/*******************************************************************************************************/
	
	if (ctrl->MRdeficit_flag == 1 && MR)
	{
		/*------------------------*/
		/* I/1. leaves */
		/*------------------------*/
		diff = diffTOTAL * ((cf->leaf_day_MR + cf->leaf_night_MR) / MR);
		
		/* calculation the difference between NSC and diff (based on available amount) */
		if (diff > 0)
		{
			/* critical NSC value: NSC pool = fixed ratio of theroretical maximum of NSC value -> but not all is available */
			NSC_crit  = (epc->prop_NSCcrit) * (cs->leafc * epc->prop_NSCvsSC);
			NSC_avail = (cs->leafc_storage + cs->leafc_transfer) - NSC_crit;
			
			if (NSC_avail < 0) NSC_avail = 0;

		
			/* calculation of deficit regarding to structural/non-structural biomass  */	 
			if (diff <= NSC_avail)
			{
				diffSC  = 0;
				diffNSC = diff; 
			}
			else
			{
				diffSC  = diff - NSC_avail;
				diffNSC = NSC_avail;
			}
			
			/* calculation of flxues from NSC pools */
			if (cs->leafc_storage + cs->leafc_transfer > 0)
			{
				cf->leafc_storage_to_MRdef		 = diffNSC * cs->leafc_storage/(cs->leafc_storage + cs->leafc_transfer);
				cf->leafc_transfer_to_MRdef	 = diffNSC * cs->leafc_transfer/(cs->leafc_storage + cs->leafc_transfer);
			}

			/* calculation of flxues from SC pools */
			cf->leafc_to_MRdef = diffSC; 
			if (cf->leafc_to_MRdef > cs->leafc)
			{
				excess             = cf->leafc_to_MRdef - cs->leafc;
				cf->leafc_to_MRdef    = cs->leafc;	
				cf->leaf_day_MR   -= excess *(cf->leaf_day_MR/(cf->leaf_day_MR +cf->leaf_night_MR));
				cf->leaf_night_MR -= excess *(cf->leaf_night_MR/(cf->leaf_day_MR +cf->leaf_night_MR));
				if (!ctrl->limitMR_flag) ctrl->limitMR_flag = 1;
			}
		}	

		/*------------------------*/
		/* I/2. froot */
		/*------------------------*/
		diff = diffTOTAL * ((cf->froot_MR) / MR);
		
		/* calculation the difference between NSC and diff (based on available amount) */
		if (diff > 0)
		{
			/* critical NSC value: NSC pool = fixed ratio of theroretical maximum of NSC value -> but not all is available */
			NSC_crit  = (epc->prop_NSCcrit) * (cs->frootc * epc->prop_NSCvsSC);
			NSC_avail = (cs->frootc_storage + cs->frootc_transfer) - NSC_crit;
			
			if (NSC_avail < 0) NSC_avail = 0;

	    	/* calculation of deficit regarding to structural/non-structural biomass  */	 
			if (diff <= NSC_avail)
			{
				diffSC  = 0;
				diffNSC = diff; 
			}
			else
			{
				diffSC  = diff - NSC_avail;
				diffNSC = NSC_avail;
			}
			
			/* calculation of flxues from NSC pools */
			if (cs->frootc_storage + cs->frootc_transfer > 0)
			{
				cf->frootc_storage_to_MRdef   = diffNSC * cs->frootc_storage/(cs->frootc_storage + cs->frootc_transfer);
				cf->frootc_transfer_to_MRdef  = diffNSC * cs->frootc_transfer/(cs->frootc_storage + cs->frootc_transfer);
			}

			/* calculation of flxues from SC pools */
			cf->frootc_to_MRdef = diffSC; 
			if (cf->frootc_to_MRdef > cs->frootc)
			{
				excess            = cf->frootc_to_MRdef - cs->frootc;
				cf->frootc_to_MRdef  = cs->frootc;	
				cf->froot_MR     -= excess;
				if (!ctrl->limitMR_flag) ctrl->limitMR_flag = 1;
			}
		}	

		/*------------------------*/
		/* I/3. softstem */
		/*------------------------*/
		diff = diffTOTAL * ((cf->softstem_MR) / MR);
		
		/* calculation the difference between NSC and diff (based on available amount) */
		if (diff > 0)
		{
			/* critical NSC value: NSC pool = fixed ratio of theroretical maximum of NSC value -> but not all is available */
			NSC_crit  = (epc->prop_NSCcrit) * (cs->softstemc * epc->prop_NSCvsSC);
			NSC_avail = (cs->softstemc_storage + cs->softstemc_transfer) - NSC_crit;
			
			if (NSC_avail < 0) NSC_avail = 0;

			/* calculation of deficit regarding to structural/non-structural biomass  */	 
			if (diff <= NSC_avail)
			{
				diffSC  = 0;
				diffNSC = diff; 
			}
			else
			{
				diffSC  = diff - NSC_avail;
				diffNSC = NSC_avail;
			}
			
			/* calculation of flxues from NSC pools */
			if (cs->softstemc_storage + cs->softstemc_transfer > 0)
			{
				cf->softstemc_storage_to_MRdef		 = diffNSC * cs->softstemc_storage/(cs->softstemc_storage + cs->softstemc_transfer);
				cf->softstemc_transfer_to_MRdef	 = diffNSC * cs->softstemc_transfer/(cs->softstemc_storage + cs->softstemc_transfer);
			}

			/* calculation of flxues from SC pools */
			cf->softstemc_to_MRdef = diffSC; 
			if (cf->softstemc_to_MRdef > cs->softstemc)
			{
				excess               = cf->softstemc_to_MRdef - cs->softstemc;
				cf->softstemc_to_MRdef  = cs->softstemc;	
				cf->softstem_MR     -= excess;
				if (!ctrl->limitMR_flag) ctrl->limitMR_flag = 1;
			}
		}	

		/*------------------------*/
		/* I/4. yield */
		/*------------------------*/
		diff = diffTOTAL * ((cf->yield_MR) / MR);
		
		/* calculation the difference between NSC and diff (based on available amount) */
		if (diff > 0)
		{
			/* critical NSC value: NSC pool = fixed ratio of theroretical maximum of NSC value -> but not all is available */
			NSC_crit  = (epc->prop_NSCcrit) * (cs->yieldc * epc->prop_NSCvsSC);
			NSC_avail = (cs->yieldc_storage + cs->yieldc_transfer) - NSC_crit;
			
			if (NSC_avail < 0) NSC_avail = 0;

			/* calculation of deficit regarding to structural/non-structural biomass  */	 
			if (diff <= NSC_avail)
			{
				diffSC  = 0;
				diffNSC = diff; 
			}
			else
			{
				diffSC  = diff - NSC_avail;
				diffNSC = NSC_avail;
			}
			
			/* calculation of flxues from NSC pools */
			if (cs->yieldc_storage + cs->yieldc_transfer > 0)
			{
				cf->yieldc_storage_to_MRdef	= diffNSC * cs->yieldc_storage/(cs->yieldc_storage + cs->yieldc_transfer);
				cf->yieldc_transfer_to_MRdef	= diffNSC * cs->yieldc_transfer/(cs->yieldc_storage + cs->yieldc_transfer);
			}

			/* calculation of flxues from SC pools */
			cf->yieldc_to_MRdef = diffSC; 
			if (cf->yieldc_to_MRdef > cs->yieldc)
			{
				excess            = cf->yieldc_to_MRdef - cs->yieldc;
				cf->yieldc_to_MRdef  = cs->yieldc;	
				cf->yield_MR     -= excess;
				if (!ctrl->limitMR_flag) ctrl->limitMR_flag = 1;
			}
		}	

		/*------------------------*/
		/* I/5. livestem */
		/*------------------------*/
		diff = diffTOTAL * ((cf->livestem_MR) / MR);
		
		/* calculation the difference between NSC and diff (based on available amount) */
		if (diff > 0)
		{
			/* critical NSC value: NSC pool = fixed ratio of theroretical maximum of NSC value -> but not all is available */
			NSC_crit  = (epc->prop_NSCcrit) * (cs->livestemc * epc->prop_NSCvsSC);
			NSC_avail = (cs->livestemc_storage + cs->livestemc_transfer) - NSC_crit;
			
			if (NSC_avail < 0) NSC_avail = 0;

		   /* calculation of deficit regarding to structural/non-structural biomass  */	 
			if (diff <= NSC_avail)
			{
				diffSC  = 0;
				diffNSC = diff; 
			}
			else
			{
				diffSC  = diff - NSC_avail;
				diffNSC = NSC_avail;
			}
			
			/* calculation of flxues from NSC pools */
			if (cs->livestemc_storage + cs->livestemc_transfer > 0)
			{
				cf->livestemc_storage_to_MRdef		 = diffNSC * cs->livestemc_storage/(cs->livestemc_storage + cs->livestemc_transfer);
				cf->livestemc_transfer_to_MRdef	 = diffNSC * cs->livestemc_transfer/(cs->livestemc_storage + cs->livestemc_transfer);
			}

			/* calculation of flxues from SC pools */
			cf->livestemc_to_MRdef = diffSC; 
			if (cf->livestemc_to_MRdef > cs->livestemc)
			{
				excess               = cf->livestemc_to_MRdef - cs->livestemc;
				cf->livestemc_to_MRdef  = cs->livestemc;	
				cf->livestem_MR     -= excess;
				if (!ctrl->limitMR_flag) ctrl->limitMR_flag = 1;
			}
		}	

		/*------------------------*/
		/* I/6. livecroot */
		/*------------------------*/
		diff = diffTOTAL * ((cf->livecroot_MR) / MR);
		
		/* calculation the difference between NSC and diff (based on available amount) */
		if (diff > 0)
		{
			/* critical NSC value: NSC pool = fixed ratio of theroretical maximum of NSC value -> but not all is available */
			NSC_crit  = (epc->prop_NSCcrit) * (cs->livecrootc * epc->prop_NSCvsSC);
			NSC_avail = (cs->livecrootc_storage + cs->livecrootc_transfer) - NSC_crit;
			
			if (NSC_avail < 0) NSC_avail = 0;

		    /* calculation of deficit regarding to structural/non-structural biomass  */	 
			if (diff <= NSC_avail)
			{
				diffSC  = 0;
				diffNSC = diff; 
			}
			else
			{
				diffSC  = diff - NSC_avail;
				diffNSC = NSC_avail;
			}
			
			/* calculation of flxues from NSC pools */
			if (cs->livecrootc_storage + cs->livecrootc_transfer > 0)
			{
				cf->livecrootc_storage_to_MRdef	 = diffNSC * cs->livecrootc_storage/(cs->livecrootc_storage + cs->livecrootc_transfer);
				cf->livecrootc_transfer_to_MRdef	 = diffNSC * cs->livecrootc_transfer/(cs->livecrootc_storage + cs->livecrootc_transfer);
			}

			/* calculation of flxues from SC pools */
			cf->livecrootc_to_MRdef = diffSC; 
			if (cf->livecrootc_to_MRdef > cs->livecrootc)
			{		
				excess                = cf->livecrootc_to_MRdef - cs->livecrootc;
				cf->livecrootc_to_MRdef  = cs->livecrootc;	
				cf->livecroot_MR     -= excess;
				if (!ctrl->limitMR_flag) ctrl->limitMR_flag = 1;
			}
		}	

	}

	
	/********************************************************************************************************/
	/* IIÍI: method 2.: woody/non-woody separated MR */
	/*******************************************************************************************************/

	if (ctrl->MRdeficit_flag == 2 && MR)
	{
		diffTOTAL_nw = diffTOTAL * (MR_nw / MR);
		diffTOTAL_w  = diffTOTAL * (MR_w  / MR);

		/*------------------------*/
		/* II/1. non-woody */
		/*------------------------*/
		/* calculation the difference between NSC and diff (based on available amount) */
		if (diffTOTAL_nw > 0)
		{
			/* critical NSC value: NSC pool = fixed ratio of theroretical maximum of NSC value -> but not all is available */
			NSC_crit  = (epc->prop_NSCcrit) * (SCnw * epc->prop_NSCvsSC);
			NSC_avail = NSCnw - NSC_crit;
			
			if (NSC_avail < 0) NSC_avail = 0;

			/* calculation of deficit regarding to structural biomass  */
			if (diffTOTAL_nw <= NSC_avail)
			{
				diffNSC = diffTOTAL_nw;
				diffSC  = 0;
			}
			else
			{
				diffNSC = NSC_avail;
				diffSC  = diffTOTAL_nw - NSC_avail;		
			}
			
			/* calculation of flxues from NSC pools */
			if (NSCnw)
			{
				cf->leafc_storage_to_MRdef		 = diffNSC * cs->leafc_storage/NSCnw;
				cf->frootc_storage_to_MRdef	 = diffNSC * cs->frootc_storage/NSCnw;
				cf->yieldc_storage_to_MRdef	 = diffNSC * cs->yieldc_storage/NSCnw;
				cf->softstemc_storage_to_MRdef	 = diffNSC * cs->softstemc_storage/NSCnw;

				cf->leafc_transfer_to_MRdef	 = diffNSC * cs->leafc_transfer/NSCnw;
				cf->frootc_transfer_to_MRdef	 = diffNSC * cs->frootc_transfer/NSCnw;
				cf->yieldc_transfer_to_MRdef	 = diffNSC * cs->yieldc_transfer/NSCnw;
				cf->softstemc_transfer_to_MRdef = diffNSC * cs->softstemc_transfer/NSCnw;	
			}
			else
			{
				if (diffNSC) 
				{
					printf("ERROR in MRdeficit calculaton in state_update.c\n");
					errorCode = 1;
				}
			}

			/* if NSC is not enough -> transfer from actual pool */
			if (diffSC > 0)
			{			
				cf->leafc_to_MRdef = diffSC * ((cf->leaf_day_MR+cf->leaf_night_MR) / MR_nw); 
				if (cf->leafc_to_MRdef > cs->leafc)
				{
					cf->leafc_to_MRdef = cs->leafc;
					excess += diffSC * ((cf->leaf_day_MR+cf->leaf_night_MR) / MR_nw) - cf->leafc_to_MRdef;
				}


				cf->frootc_to_MRdef = diffSC * (cf->froot_MR / MR_nw); 
				if (cf->frootc_to_MRdef > cs->frootc)
				{
					cf->frootc_to_MRdef = cs->frootc;
					excess += diffSC * (cf->froot_MR / MR_nw) - cf->frootc_to_MRdef;
				}


				cf->yieldc_to_MRdef = diffSC * (cf->yield_MR / MR_nw); 
				if (cf->yieldc_to_MRdef > cs->yieldc)
				{
					cf->yieldc_to_MRdef = cs->yieldc;
					excess += diffSC * (cf->yield_MR / MR_nw) - cf->yieldc_to_MRdef;
				}

				cf->softstemc_to_MRdef = diffSC * (cf->softstem_MR / MR_nw); 
				if (cf->softstemc_to_MRdef > cs->softstemc)
				{
					cf->softstemc_to_MRdef = cs->softstemc;
					excess += diffSC * (cf->softstem_MR / MR_nw) - cf->softstemc_to_MRdef;
				}


				/* if maintresp of non-woody biomass can not be covered from non-woody biomass -> decreasing of maintanance flux */
				if (excess)
				{
					if (!ctrl->limitMR_flag) ctrl->limitMR_flag = 1;
					cs->MRlimitSUM    += excess;
					reduction         += excess * (cf->leaf_day_MR   / MR_nw);
					cf->leaf_day_MR   -= excess * (cf->leaf_day_MR   / MR_nw);
					reduction         += excess * (cf->leaf_night_MR / MR_nw);
 					cf->leaf_night_MR -= excess * (cf->leaf_night_MR / MR_nw);
					reduction         += excess * (cf->froot_MR      / MR_nw);
					cf->froot_MR      -= excess * (cf->froot_MR      / MR_nw);
					reduction         += excess * (cf->softstem_MR   / MR_nw);
					cf->softstem_MR   -= excess * (cf->softstem_MR   / MR_nw);
					reduction         += excess * (cf->yield_MR      / MR_nw);
					cf->yield_MR      -= excess * (cf->yield_MR      / MR_nw);
					excess            -= reduction;
					if (fabs(excess) > CRIT_PREC)
					{
						printf("ERROR in MRdeficit calculaton in state_update.c\n");
						errorCode = 1;
					}
				}
			}
		}

		/*------------------------*/
		/* II/1. woody */
		/*------------------------*/
	
		if (MR_w)
		{
			/* 2.1. calculation the difference between NSC and diff (based on available amount) */
			if (diffTOTAL_w > 0)
			{
				/* critical NSC value: NSC pool = fixed ratio of theroretical maximum of NSC value -> but not all is available */
				NSC_crit  = (epc->prop_NSCcrit) * (SCw * epc->prop_NSCvsSC);
				NSC_avail = NSCw - NSC_crit;

				if (NSC_avail < 0) NSC_avail = 0;
		
				/* calculation of deficit regarding to structural biomass  */
 				if (diffTOTAL_w <= NSC_avail)
				{
					diffNSC = diffTOTAL_w;
					diffSC  = 0;
				}
				else
				{
					diffNSC = NSC_avail;
					diffSC  = diffTOTAL_w - NSC_avail;		
				}

				/* 2.2. calculation of flxues from nsc pools */
				if (NSCw)
				{
					cf->livestemc_storage_to_MRdef	 = diffNSC * cs->livestemc_storage/NSCw;
					cf->livecrootc_storage_to_MRdef = diffNSC * cs->livecrootc_storage/NSCw;
					cf->deadstemc_storage_to_MRdef	 = diffNSC * cs->deadstemc_storage/NSCw;
					cf->deadcrootc_storage_to_MRdef = diffNSC * cs->deadcrootc_storage/NSCw;

					cf->livestemc_transfer_to_MRdef  = diffNSC * cs->livestemc_transfer/NSCw;
					cf->livecrootc_transfer_to_MRdef = diffNSC * cs->livecrootc_transfer/NSCw;
					cf->deadstemc_transfer_to_MRdef  = diffNSC * cs->deadstemc_transfer/NSCw;
					cf->deadcrootc_transfer_to_MRdef = diffNSC * cs->deadcrootc_transfer/NSCw;			
				}
				else
				{
					if (diffNSC) 
					{
						printf("ERROR in MRdeficit calculaton in state_update.c\n");
						errorCode = 1;
					}
				}

				/* 1.5. if NSC is not enough -> transfer from actual pool */
				if (diffSC > 0)
				{
					cf->livestemc_to_MRdef = diffSC * (cf->livestem_MR / MR_w); 	
					if (cs->livestemc < cf->livestemc_to_MRdef)
					{
						cf->livestemc_to_MRdef = cs->livestemc;
						excess += diffSC * (cf->livestem_MR / MR_w) - cf->livestemc_to_MRdef;
					}

					cf->livecrootc_to_MRdef = diffSC * (cf->livecroot_MR / MR_w);
					if (cs->livecrootc < cf->livecrootc_to_MRdef)
					{
						cf->livecrootc_to_MRdef = cs->livecrootc;
						excess += diffSC * (cf->livecroot_MR / MR_w) - cf->livecrootc_to_MRdef;
					}

					/* if maintresp of woody biomass can not be covered from non-woody biomass -> decreasing of maintanance flux */
					if (excess)
					{
						cs->MRlimitSUM   += excess;
						cf->livestem_MR  -= excess * cf->livestem_MR  / (MR_w);
 						cf->livecroot_MR -= excess * cf->livecroot_MR / (MR_w);
						excess             = 0;
						if (!ctrl->limitMR_flag) ctrl->limitMR_flag = 1;
					}

				}
			}
		}
	}


	if (ctrl->MRdeficit_flag != 3)
	{
		/* calculation of N-fluxes */
		if (epc->leaf_cn)     
		{
			nf->leafn_storage_to_MRdef		  = cf->leafc_storage_to_MRdef / epc->leaf_cn;
			nf->leafn_transfer_to_MRdef      = cf->leafc_transfer_to_MRdef / epc->leaf_cn;
			nf->leafn_to_MRdef               = (cf->leafc_to_MRdef)  / epc->leaf_cn;
		}
		if (epc->froot_cn)    
		{
			nf->frootn_storage_to_MRdef	   = cf->frootc_storage_to_MRdef / epc->froot_cn;
			nf->frootn_transfer_to_MRdef      = cf->frootc_transfer_to_MRdef / epc->froot_cn;
			nf->frootn_to_MRdef               = cf->frootc_to_MRdef / epc->froot_cn;
		}
		if (epc->yield_cn)    
		{
			nf->yieldn_storage_to_MRdef		= cf->yieldc_storage_to_MRdef / epc->yield_cn;
			nf->yieldn_transfer_to_MRdef		= cf->yieldc_transfer_to_MRdef / epc->yield_cn;
			nf->yieldn_to_MRdef                = cf->yieldc_to_MRdef / epc->yield_cn;
		}
					
		if (epc->softstem_cn) 
		{
			nf->softstemn_storage_to_MRdef      = cf->softstemc_storage_to_MRdef / epc->softstem_cn;
			nf->softstemn_transfer_to_MRdef     = cf->softstemc_transfer_to_MRdef / epc->softstem_cn;
			nf->softstemn_to_MRdef              = cf->softstemc_to_MRdef / epc->softstem_cn;
		}

				
		if (epc->livewood_cn)
		{
			nf->livestemn_storage_to_MRdef   = cf->livestemc_storage_to_MRdef / epc->livewood_cn;
			nf->livecrootn_storage_to_MRdef  = cf->livecrootc_storage_to_MRdef / epc->livewood_cn;
			nf->livestemn_transfer_to_MRdef  = cf->livestemc_transfer_to_MRdef / epc->livewood_cn;
			nf->livecrootn_transfer_to_MRdef = cf->livecrootc_transfer_to_MRdef / epc->livewood_cn;
			nf->livestemn_to_MRdef           = cf->livestemc_to_MRdef / epc->livewood_cn; 
			nf->livecrootn_to_MRdef          = cf->livecrootc_to_MRdef / epc->livewood_cn;		
			nf->livestemn_to_MRdef           = cf->livestemc_to_MRdef / epc->livewood_cn;
			nf->livecrootn_to_MRdef          = cf->livecrootc_to_MRdef / epc->livewood_cn;
		}

		if (epc->deadwood_cn) 
		{
			nf->deadstemn_storage_to_MRdef   = cf->deadstemc_storage_to_MRdef / epc->deadwood_cn;
			nf->deadcrootn_storage_to_MRdef  = cf->deadcrootc_storage_to_MRdef / epc->deadwood_cn;
			nf->deadstemn_transfer_to_MRdef  = cf->deadstemc_transfer_to_MRdef / epc->deadwood_cn;
			nf->deadcrootn_transfer_to_MRdef = cf->deadcrootc_transfer_to_MRdef / epc->deadwood_cn;
		}
	
		/* summarize fluxes */		
		cf->NSCnw_to_MRdef = cf->leafc_storage_to_MRdef  + cf->frootc_storage_to_MRdef  + cf->yieldc_storage_to_MRdef  + cf->softstemc_storage_to_MRdef +
						  cf->leafc_transfer_to_MRdef + cf->frootc_transfer_to_MRdef + cf->yieldc_transfer_to_MRdef + cf->softstemc_transfer_to_MRdef;

		cf->NSCw_to_MRdef  = cf->livestemc_storage_to_MRdef  + cf->livecrootc_storage_to_MRdef  + cf->deadstemc_storage_to_MRdef  + cf->deadcrootc_storage_to_MRdef +
						  cf->livestemc_transfer_to_MRdef + cf->livecrootc_transfer_to_MRdef + cf->deadstemc_transfer_to_MRdef + cf->deadcrootc_transfer_to_MRdef;
	
		cf->SCnw_to_MRdef  = cf->leafc_to_MRdef + cf->frootc_to_MRdef + cf->yieldc_to_MRdef + cf->softstemc_to_MRdef;

		cf->SCw_to_MRdef   = cf->livestemc_to_MRdef + cf->livecrootc_to_MRdef;

		nf->NSNnw_to_MRdef = nf->leafn_storage_to_MRdef  + nf->frootn_storage_to_MRdef  + nf->yieldn_storage_to_MRdef  + nf->softstemn_storage_to_MRdef +
						  nf->leafn_transfer_to_MRdef + nf->frootn_transfer_to_MRdef + nf->yieldn_transfer_to_MRdef + nf->softstemn_transfer_to_MRdef;
	
		nf->NSNw_to_MRdef  = nf->livestemn_storage_to_MRdef  + nf->livecrootn_storage_to_MRdef  + nf->deadstemn_storage_to_MRdef  + nf->deadcrootn_storage_to_MRdef +
						  nf->livestemn_transfer_to_MRdef + nf->livecrootn_transfer_to_MRdef + nf->deadstemn_transfer_to_MRdef + nf->deadcrootn_transfer_to_MRdef;
			
		nf->SNnw_to_MRdef  = nf->leafn_to_MRdef + nf->frootn_to_MRdef + nf->yieldn_to_MRdef + nf->softstemn_to_MRdef;	

		nf->SNw_to_MRdef   = nf->livestemn_to_MRdef + nf->livecrootn_to_MRdef;

	

		/* state update of storage and transfer pools */
		cs->leafc_storage					-= cf->leafc_storage_to_MRdef;
		cs->frootc_storage					-= cf->frootc_storage_to_MRdef;
		cs->yieldc_storage					-= cf->yieldc_storage_to_MRdef;
		cs->softstemc_storage				-= cf->softstemc_storage_to_MRdef;
		
		cs->leafc_transfer					-= cf->leafc_transfer_to_MRdef;
		cs->frootc_transfer					-= cf->frootc_transfer_to_MRdef;
		cs->yieldc_transfer					-= cf->yieldc_transfer_to_MRdef;
		cs->softstemc_transfer				-= cf->softstemc_transfer_to_MRdef;
		
		ns->leafn_storage					-= nf->leafn_storage_to_MRdef;
		ns->frootn_storage					-= nf->frootn_storage_to_MRdef;
		ns->yieldn_storage					-= nf->yieldn_storage_to_MRdef;
		ns->softstemn_storage				-= nf->softstemn_storage_to_MRdef;
		
		ns->leafn_transfer					-= nf->leafn_transfer_to_MRdef;
		ns->frootn_transfer					-= nf->frootn_transfer_to_MRdef;
		ns->yieldn_transfer					-= nf->yieldn_transfer_to_MRdef;
		ns->softstemn_transfer				-= nf->softstemn_transfer_to_MRdef;
	
		cs->livestemc_storage				-= cf->livestemc_storage_to_MRdef;
		cs->livecrootc_storage				-= cf->livecrootc_storage_to_MRdef;
		cs->deadstemc_storage				-= cf->deadstemc_storage_to_MRdef;
		cs->deadcrootc_storage				-= cf->deadcrootc_storage_to_MRdef;

		cs->livestemc_transfer				-= cf->livestemc_transfer_to_MRdef;
		cs->livecrootc_transfer				-= cf->livecrootc_transfer_to_MRdef;
		cs->deadstemc_transfer				-= cf->deadstemc_transfer_to_MRdef;
		cs->deadcrootc_transfer				-= cf->deadcrootc_transfer_to_MRdef;
			
		ns->livestemn_storage				-= nf->livestemn_storage_to_MRdef;
		ns->livecrootn_storage				-= nf->livecrootn_storage_to_MRdef;
		ns->deadstemn_storage				-= nf->deadstemn_storage_to_MRdef;
		ns->deadcrootn_storage				-= nf->deadcrootn_storage_to_MRdef;

		ns->livestemn_transfer				-= nf->livestemn_transfer_to_MRdef;
		ns->livecrootn_transfer				-= nf->livecrootn_transfer_to_MRdef;
		ns->deadstemn_transfer				-= nf->deadstemn_transfer_to_MRdef;
		ns->deadcrootn_transfer				-= nf->deadcrootn_transfer_to_MRdef;
	
	
		/* state update of actual pools */
		cs->leafc					-= cf->leafc_to_MRdef;
		cs->frootc					-= cf->frootc_to_MRdef;
		cs->yieldc					-= cf->yieldc_to_MRdef;
		cs->softstemc				-= cf->softstemc_to_MRdef;

		ns->leafn					-= nf->leafn_to_MRdef;
		ns->frootn					-= nf->frootn_to_MRdef;
		ns->yieldn					-= nf->yieldn_to_MRdef;
		ns->softstemn				-= nf->softstemn_to_MRdef;
			
		cs->livestemc				-= cf->livestemc_to_MRdef;
		cs->livecrootc				-= cf->livecrootc_to_MRdef;

		ns->livestemn				-= nf->livestemn_to_MRdef;
		ns->livecrootn				-= nf->livecrootn_to_MRdef;

		ns->retransn                += nf->NSNnw_to_MRdef + nf->NSNw_to_MRdef + nf->SNnw_to_MRdef + nf->SNw_to_MRdef;
				
	}

	MR_nw = cf->leaf_day_MR + cf->leaf_night_MR + cf->froot_MR + cf->yield_MR + cf->softstem_MR;
	MR_w  = cf->livestem_MR + cf->livecroot_MR;
	MR    = MR_nw + MR_w;

	/* control */			
	diff=(MR - (cs->cpool+cf->NSCnw_to_MRdef+cf->SCnw_to_MRdef+cf->NSCw_to_MRdef+cf->SCw_to_MRdef));
	if (ctrl->MRdeficit_flag != 3 && diff > CRIT_PREC && diff+cs->cpool > CRIT_PREC)
	{
  		printf("ERROR in maintanance respiration deficit calculation in state_update.c\n");
		errorCode=1;
	}

	/* state update of cpool and npool (the released N turns into retranslocated N pool) */
	if (ns->npool) CpoolNpool_ratio = cs->cpool / ns->npool;
	
	cf->cpool_to_MRdef  = (MR - cf->NSCnw_to_MRdef - cf->SCnw_to_MRdef - cf->NSCw_to_MRdef - cf->SCw_to_MRdef);

	
	cs->cpool	   -= cf->cpool_to_MRdef;

	if (CpoolNpool_ratio && ns->npool > 0)
	{
		npool_to_MRdef = cf->cpool_to_MRdef / CpoolNpool_ratio;
		if (ns->npool < npool_to_MRdef) npool_to_MRdef = ns->npool;
		ns->npool     -= npool_to_MRdef;
		ns->retransn += npool_to_MRdef;
	}

		
	/* state update MR sink pools */	
	cs->MRleaf_snk			+= cf->leaf_day_MR;
	cs->MRleaf_snk			+= cf->leaf_night_MR;
	cs->MRfroot_snk			+= cf->froot_MR;
	cs->MRyield_snk			+= cf->yield_MR;
	cs->MRsoftstem_snk	    += cf->softstem_MR;
	cs->MRlivestem_snk		+= cf->livestem_MR;
	cs->MRlivecroot_snk		+= cf->livecroot_MR;

		
	cs->NSCnw = cs->leafc_storage      +  cs->frootc_storage     + cs->yieldc_storage     + cs->softstemc_storage +
					cs->leafc_transfer     +  cs->frootc_transfer    + cs->yieldc_transfer    + cs->softstemc_transfer;

	cs->NSCw = cs->livestemc_storage  + cs->livecrootc_storage  + cs->deadstemc_storage  + cs->deadcrootc_storage +
				cs->livestemc_transfer + cs->livecrootc_transfer + cs->deadstemc_transfer + cs->deadcrootc_transfer;

	cs->SCnw = cs->leafc     +  cs->frootc    + cs->yieldc     + cs->softstemc;

	cs->SCw = cs->livestemc  + cs->livecrootc  + cs->deadstemc  + cs->deadcrootc;


	return (errorCode);
}			

