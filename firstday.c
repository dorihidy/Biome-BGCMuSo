/* 
firstday.c
Initializes the state variables for the first day of a simulation that is not using a restart file.

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGCMuSo v7.0.
Original code: Copyright 2000, Peter E. Thornton
Numerical Terradynamic Simulation Group, The University of Montana, USA
Modified code: Copyright 2025, D. Hidy [dori.hidy@gmail.com]
Hungarian Academy of Sciences, Hungary
See the website of Biome-BGCMuSo at http://nimbus.elte.hu/bbgc/ for documentation, model executable and example input files.
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include "ini.h"
#include "bgc_struct.h"     /* structure definitions */
#include "bgc_func.h"       /* function prototypes */
#include "bgc_constants.h"


int firstday(const control_struct* ctrl, const epconst_struct* epc, const planting_struct* PLT, 
	         soilprop_struct* sprop, siteconst_struct* sitec, cinit_struct* cinit, phenology_struct* phen, epvar_struct* epv, soilInfo_struct* soilInfo,
	         wstate_struct* ws, cstate_struct* cs, nstate_struct* ns, summary_struct* summary, psn_struct* psn_sun, psn_struct* psn_shade)
{
	int errorCode=0;
	int layer, day, pp, dm;
	double prop_transfer, transfer, prop_litfall;
	double max_deadstemc, max_deadcrootc;
	int N_DECOMPLAYERS = 3;
	
	/* *****************************************************************************- */
	/* 1. Initialize ecophysiological variables */

	epv->DSR = 0.0;
	epv->cumWS_anoxic = 0;
	epv->cumWS_drought = 0;
	epv->cumWS = 0.0;
	epv->cumNS = 0.0;
	epv->WSlenght = 0.0;
	epv->transfer_ratio = 0.0;
	epv->leafday = -1;
	epv->leafday_lastmort = -1;

	epv->n_rootlayers = 0;
	epv->germ_layer = 1;

	epv->germDepth = 0.05;
    epv->projLAI = 0;
	epv->projLAI_STDB = 0;
    epv->allLAI = 0;
	epv->SLA_avg = 0;
    epv->plaisun = 0;
    epv->plaishade = 0;
    epv->projSLA_sun = 0;
    epv->projSLA_shade = 0;
	epv->plantHeight = 0;
	epv->n_actphen = 0;
	epv->flowHS_mort = 0;
	epv->assim_sun = 0;
	epv->assim_shade = 0;
	epv->dlmr_area_sun = 0;
	epv->dlmr_area_shade = 0;

	epv->plantCalloc = 0; 
	epv->plantNalloc = 0;
    epv->cumCalloc_plant = 0; 
	epv->cumNalloc_plant = 0;
	epv->excess_c = 0;
	epv->pnow = 0;
	epv->NDVI = 0;

	epv->rootlength = 0;
	epv->rs_decomp_avg = 0;
	epv->m_Tmin = 0;
	epv->m_co2 = 0;
	epv->stomaCONDUCT_max = 0;
	epv->albedo_LAI = 0;
	epv->assim_Tcoeff = 1;
	epv->assim_WScoeff = 1;
	epv->SCpercent = 0;
	epv->SC_EVPred = 1;
	epv->litr1_CNratio = 0;
	epv->litr2_CNratio = 0;
	epv->litr3_CNratio = 0;
	epv->litr4_CNratio = 0;
					
	psn_sun->A      = 0;
	psn_sun->Ci	    = 0;
	psn_sun->O2	    = 0;
	psn_sun->Ca	    = 0;
	psn_sun->gamma	= 0;
	psn_sun->Kc	    = 0;
	psn_sun->Ko	    = 0;
	psn_sun->Vmax	= 0;
	psn_sun->Jmax	= 0;
	psn_sun->J	    = 0;
	psn_sun->Av	    = 0;
	psn_sun->Aj	    = 0;
	
	psn_shade->A       = 0;
	psn_shade->Ci	    = 0;
	psn_shade->O2	    = 0;
	psn_shade->Ca	    = 0;
	psn_shade->gamma	= 0;
	psn_shade->Kc	    = 0;
	psn_shade->Ko	    = 0;
	psn_shade->Vmax	    = 0;
	psn_shade->Jmax  	= 0;
	psn_shade->J	    = 0;
	psn_shade->Av	    = 0;
	psn_shade->Aj	    = 0;

	
	summary->leafDM = 0;
	summary->leaflitrDM = 0;
	summary->frootDM = 0;
	summary->yieldDM = 0;
	summary->softstemDM = 0;
	summary->livewoodDM = 0;
	summary->deadwoodDM = 0;
	summary->yieldDM_HRV = 0;
	summary->vegC = 0;
	summary->vegN = 0;
	summary->orgN = 0;
	summary->LDaboveC_nw = 0;
	summary->LDaboveC_w = 0;
	summary->LDaboveCwithNSC_nw = 0;
	summary->LDaboveCwithNSC_w = 0;
	summary->LaboveC_nw = 0;
	summary->LaboveC_w = 0;
	summary->LaboveCwithNSC_nw = 0;
	summary->LaboveCwithNSC_w = 0;
	summary->DaboveC_nw = 0;
	summary->DaboveC_w = 0;
	summary->DaboveCwithNSC_nw = 0;
	summary->DaboveCwithNSC_w = 0;
	summary->LDbelowC_nw = 0;
	summary->LDbelowC_w = 0;
	summary->LDbelowCwithNSC_nw = 0;
	summary->LDbelowCwithNSC_w = 0;
	summary->LbelowC_nw = 0;
	summary->LbelowC_w = 0;
	summary->LbelowCwithNSC_nw = 0;
	summary->LbelowCwithNSC_w = 0;
	summary->DbelowC_nw = 0;
	summary->DbelowC_w = 0;
	summary->DbelowCwithNSC_nw = 0;
	summary->DbelowCwithNSC_w = 0;
	summary->livingSC = 0;
	summary->livingNSC = 0;
	summary->livingBIOMabove = 0;
	summary->livingBIOMbelow = 0;
	summary->BIOMaboveSUM = 0;
	summary->BIOMbelowSUM = 0;
	
	summary->litrCwdC_total = 0;
	summary->litrCwdN_total = 0;
	summary->litrN_total = 0;
	summary->litrC_total = 0;
	summary->soilC_total = 0;
	summary->soilN_total = 0;
	summary->sminN_total = 0;
	summary->sminNavail_total = 0;
	summary->sminNavail_RZmax = 0;
	summary->sminN_RZmax = 0;
	summary->NO3_RZmax = 0;
	summary->NH4_RZmax = 0;
	summary->Wleach_RZmax = 0;
	summary->DOCleach_RZmax = 0;
	summary->DONleach_RZmax = 0;
	summary->sminNleach_RZmax = 0;
	summary->soilC_RZmax = 0;
	summary->soilN_RZmax = 0;
	summary->litrC_RZmax = 0;
	summary->litrN_RZmax = 0;
	summary->totalC = 0;
	summary->totalN = 0;

	
	summary->leafc_LandD = 0;
	summary->frootc_LandD = 0;
	summary->yield_LandD = 0;
	summary->softstemc_LandD = 0;
	summary->CNlitr_total = 0;
	summary->NH4_unsat = 0;
	summary->NO3_unsat = 0;
	summary->orgN_unsat = 0;



	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		/* possibility to calculate dynamic changing Ksat 		*/
		if (sprop->hydrCONTUCTsatMES_cmPERday[layer] == (double)DATA_GAP)
		{
			if (ctrl->Ksat_flag == 1)
				sprop->hydrCONDUCTsat[layer] = 7.05556 * 1e-6 * pow(10, (-0.6 + 0.0126 * sprop->sand[layer] - 0.0064 * sprop->clay[layer]));
			else
				sprop->hydrCONDUCTsat[layer] = (50 * exp(-0.075 * sprop->clay[layer]) + 200 * exp(-0.075 * (100 - sprop->sand[layer]))) / 100 / nSEC_IN_DAY;
		}

		sprop->hydrDIFFUSsat[layer] = (sprop->soilB[layer] * sprop->hydrCONDUCTsat[layer] * (-100 * sprop->PSIsat[layer])) / sprop->VWCsat[layer];
		sprop->hydrCONDUCTfc[layer] = sprop->hydrCONDUCTsat[layer] * pow(sprop->VWCfc[layer] / sprop->VWCsat[layer], 2 * sprop->soilB[layer] + 3);
		sprop->hydrDIFFUSfc[layer] = (((sprop->soilB[layer] * sprop->hydrCONDUCTsat[layer] * (-100 * sprop->PSIsat[layer]))) / sprop->VWCsat[layer])
			* pow(sprop->VWCfc[layer] / sprop->VWCsat[layer], sprop->soilB[layer] + 2);

		summary->NH4_ppm[layer] = 0;
		summary->NO3_ppm[layer] = 0;
		summary->orgN_ppm[layer] = 0;
		summary->sminNavail[layer] = 0;
		summary->SOCpercent[layer] = 0;


	}

	for (layer = 0; layer < N_SOILvirtLAYERS; layer++)
	{
		sprop->Ztop[layer] = 0;
		sprop->Zbot[layer] = 0;
		sprop->virtLayer[layer] = 0;
	}

	sprop->ratioNORMcf = 0;
	sprop->ratioCAPILcf = 0;
	sprop->ratioNORMgw = 0;
	sprop->ratioCAPILgw = 0;

	sprop->balance_NORMcf = 0;
	sprop->balance_CAPILcf = 0;
	sprop->balance_NORMgw = 0;
	sprop->balance_CAPILgw = 0;

	sprop->soilw_NORMcf_pre = 0;
	sprop->soilw_CAPILcf_pre = 0;
	sprop->soilw_NORMgw_pre = 0;
	sprop->soilw_CAPILgw_pre = 0;
	sprop->soilw_SATgw_pre = 0;

	/* the role of the layer in the decomposition */
	
	if (sprop->PROPlayerDC_mes[0] == DATA_GAP)
	{ 
		for (layer = 0; layer < N_DECOMPLAYERS; layer++) sprop->PROPlayerDC[layer] = sitec->soillayer_thickness[layer] / sitec->soillayer_depth[N_DECOMPLAYERS - 1];
		for (layer = N_DECOMPLAYERS; layer < N_SOILLAYERS; layer++) sprop->PROPlayerDC[layer] = 0;
	}
	else
		for (layer = 0; layer < N_SOILLAYERS; layer++) sprop->PROPlayerDC[layer] = sprop->PROPlayerDC_mes[layer];

	if (ctrl->spinup)
	{
		sprop->GWD_pre = DATA_GAP;
		for (dm = 0; dm < N_DISSOLVMATER; dm++)
		{
			soilInfo->content_NORMcf[dm] = 0;
			soilInfo->content_CAPILcf[dm] = 0;
			soilInfo->content_NORMgw[dm] = 0;
			soilInfo->content_CAPILgw[dm] = 0;
			soilInfo->content_SATgw[dm] = 0;
		}
	}
		

	/* initalize the number of the soil layers in which root can be found. It determines the rootzone depth (only on first day) */
	
	if (calc_nrootlayers(0, epc->rootzoneDepth_max, cs->frootc, sitec, epv))
	{
		if (!errorCode) 
		{
			printf("\n");
			printf("ERROR in calc_nrootlayers.c for multilayer_rootDepth.c\n");
		}
		errorCode=1;
	}

	/* initialize multilayer variables (first approximation: field cap.) and multipliers for stomatal limitation calculation */
	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		epv->m_WS_layer[layer]  = 1;
		epv->ws_nitrif[layer] = 0;
		epv->ps_nitrif[layer] = 0;
		epv->ts_nitrif[layer] = 0;
		epv->ws_decomp[layer] = 0;
		epv->rs_decomp[layer] = 0;
		epv->ts_decomp[layer] = 0;
		epv->IMMOBratio[layer] = 0;
	    epv->rootlengthProp[layer]     = 0;
		epv->rootlengthLandD_prop[layer]= 0;
		ns->sminN[layer]              = ns->NH4[layer] + ns->NO3[layer];
		cs->soil1c_total += cs->soil1c[layer];
		cs->soil2c_total += cs->soil2c[layer];
		cs->soil3c_total += cs->soil3c[layer];
		cs->soil4c_total += cs->soil4c[layer];
	}



    /* evergreen biome: root available also in the first day */
	epv->rootDepth = 0;
	if (epc->evergreen) 
	{
		epv->rootlengthProp[0] = 1;
		epv->n_rootlayers = 1;
		if (epc->woody) 
			epv->rootDepth = epc->rootzoneDepth_max;
		else
			epv->rootDepth = CRIT_PREC;
	}
	/* initialize genetical senescence variables */

	epv->thermal_time = 0;

	epv->m_ppfd_sun = 1;
	epv->m_ppfd_shade = 1;
	epv->m_vpd = 1;
	epv->m_final_sun = 1;
	epv->m_final_shade = 1;
	epv->m_WSlenght = 1;
	epv->m_extremT = 1;
	epv->gcorr = 0;
	epv->gl_bl = 0;
	epv->gl_c = 0;
	epv->gl_s_sun = 0;
	epv->gl_s_shade = 0;
	epv->gl_e_wv = 0;
	epv->gl_sh = 0;
	epv->gc_e_wv = 0;
	epv->gc_sh = 0;
	epv->gl_t_wv_sun = 0;
	epv->gl_t_wv_shade = 0;
	epv->gl_t_wv_sunPOT = 0;
	epv->gl_t_wv_shadePOT = 0;
	epv->WPM = 0;
	epv->FM = 0;
	epv->plantNdemand = 0;
	epv->ES = 0;
	epv->EV = 0;
	epv->Tday = 0;

	for (day = 0; day < 2*nDAYS_OF_YEAR; day++)
	{
		epv->thermal_timeARRAY[day]     = 0;
		epv->cpool_to_leafcARRAY[day] = 0;
		epv->npool_to_leafnARRAY[day] = 0;
		epv->gpSNSC_phenARRAY[day] = 0;
	}

	epv->m_WS	            = 1;
	epv->m_WSanoxic         = 1;
	epv->m_WSdrought        = 1;
	epv->SMSI               = 0;
	epv->flower_date        = 0;
	epv->winterEnd_date     = 0;

	phen->GDD_emergSTART = 0;
	phen->GDD_emergEND   = 0;
	phen->GDD_limit      = 0;
	phen->yday_total     = 0;
	phen->onday          = -1;
	phen->offday         = -1;
	phen->planttype      = 0;



	/* initialize sum of leafC content and starting date of a given phenophase */
	for (pp = 0; pp < N_PHENPHASES; pp++) 
	{
		cs->leafcSUM_phenphase[pp] = 0;
		epv->phenphase_date[pp]    =-1;
		epv->rootDepth_phen[pp]    = 0;
	}


	/* *****************************************************************************- */
	/* 3. Initialize the C and N  variables, but NOT in transient phase and NOT in normal phase without reading restart file */


	if (ctrl->read_restart == 0 && ctrl->spinup != 2)
	{
		/* all types can use annmax leafc and frootc */
		epv->annmax_leafc = 0.0;
		epv->annmax_frootc = 0.0;
		epv->annmax_yieldc = 0.0;
		epv->annmax_softstemc = 0.0;
		epv->annmax_livestemc = 0.0;
		epv->annmax_livecrootc = 0.0;

		/* initialize the c and n storage state variables, but NOT in transient phase */
		cs->leafc_storage = 0.0;
		cs->frootc_storage = 0.0;
		cs->yieldc_storage = 0.0;
		cs->softstemc_storage = 0.0;
		cs->livestemc_storage = 0.0;
		cs->deadstemc_storage = 0.0;
		cs->livecrootc_storage = 0.0;
		cs->deadcrootc_storage = 0.0;
		cs->gresp_storage = 0.0;
		cs->cpool = 0.0;
		
		/* initalization of above- and belowground litter */
		cs->litrCabove[0] = cs->litr1c[0] + cs->litr2c[0] + cs->litr3c[0] + cs->litr4c[0];
		cs->litrCbelow[0] = 0;
		cs->cwdCabove[0] = cs->cwdc[0];
		cs->cwdCbelow[0] = 0;

		for (layer = 1; layer < N_SOILLAYERS; layer++) 
		{ 
			cs->litrCabove[layer] = 0;
			cs->litrCbelow[layer] = cs->litr1c[layer] + cs->litr2c[layer] + cs->litr3c[layer] + cs->litr4c[layer]; 
			cs->cwdCabove[layer] = 0;
			cs->cwdCbelow[layer] = cs->cwdc[layer];
		}

		
		ns->leafn_storage = 0.0;
		ns->frootn_storage = 0.0;
		ns->yieldn_storage = 0.0;
		ns->softstemn_storage = 0.0;
		ns->livestemn_storage = 0.0;
		ns->deadstemn_storage = 0.0;
		ns->livecrootn_storage = 0.0;
		ns->deadcrootn_storage = 0.0;
		ns->npool = 0.0;
	
	
		/* establish the initial partitioning between displayed growth and growth ready for transfer */
		/* calculate initial leaf and froot nitrogen pools from carbon pools and user-specified initial C:N for each component */

		if (epc->leaf_cn)
		{
			cs->leafc_transfer = cinit->max_leafc * epc->nonwoody_turnover;
			cs->leafc = cinit->max_leafc - cs->leafc_transfer;
			ns->leafn_transfer = cs->leafc_transfer / epc->leaf_cn;
			ns->leafn = cs->leafc / epc->leaf_cn;
		}

		if (epc->froot_cn)
		{
			cinit->max_frootc = cinit->max_leafc * (epc->alloc_frootc[0] / epc->alloc_leafc[0]);
			cs->frootc_transfer = cinit->max_frootc * epc->nonwoody_turnover;
			cs->frootc = cinit->max_frootc - cs->frootc_transfer;
			ns->frootn_transfer = cs->frootc_transfer / epc->froot_cn;
			ns->frootn = cs->frootc / epc->froot_cn;
		}

		if (epc->yield_cn)
		{
			cinit->max_yieldc = cinit->max_leafc * (epc->alloc_yield[0] / epc->alloc_leafc[0]);
			cs->yieldc_transfer = cinit->max_yieldc * epc->nonwoody_turnover;
			cs->yieldc = cinit->max_yieldc - cs->yieldc_transfer;
			ns->yieldn_transfer = cs->yieldc_transfer / epc->yield_cn;
			ns->yieldn = cs->yieldc / epc->yield_cn;
		}

		if (epc->softstem_cn)
		{
			cinit->max_softstemc = cinit->max_leafc * (epc->alloc_softstemc[0] / epc->alloc_leafc[0]);
			cs->softstemc_transfer = cinit->max_softstemc * epc->nonwoody_turnover;
			cs->softstemc = cinit->max_softstemc - cs->softstemc_transfer;
			ns->softstemn_transfer = cs->softstemc_transfer / epc->softstem_cn;
			ns->softstemn = cs->softstemc / epc->softstem_cn;
		}


		if (epc->livewood_cn)
		{
			cinit->max_livestemc = cinit->max_leafc * (epc->alloc_livestemc[0] / epc->alloc_leafc[0]);
			cs->livestemc_transfer = cinit->max_livestemc * epc->woody_turnover;
			cs->livestemc = cinit->max_livestemc - cs->livestemc_transfer;
			ns->livestemn_transfer = cs->livestemc_transfer / epc->livewood_cn;
			ns->livestemn = cs->livestemc / epc->livewood_cn;

			cinit->max_livecrootc = cinit->max_leafc * (epc->alloc_livecrootc[0] / epc->alloc_leafc[0]);
			cs->livecrootc_transfer = cinit->max_livecrootc * epc->woody_turnover;
			cs->livecrootc = cinit->max_livecrootc - cs->livecrootc_transfer;
			ns->livecrootn_transfer = cs->livecrootc_transfer / epc->livewood_cn;
			ns->livecrootn = cs->livecrootc / epc->livewood_cn;
		}

		if (epc->deadwood_cn)
		{
			max_deadstemc = cinit->max_leafc * (epc->alloc_deadstemc[0] / epc->alloc_leafc[0]);
			cs->deadstemc_transfer = max_deadstemc * epc->woody_turnover;
			cs->deadstemc = max_deadstemc - cs->deadstemc_transfer;
			ns->deadstemn_transfer = cs->deadstemc_transfer / epc->deadwood_cn;
			ns->deadstemn = cs->deadstemc / epc->deadwood_cn;

			max_deadcrootc = cinit->max_leafc * (epc->alloc_deadcrootc[0] / epc->alloc_leafc[0]);
			cs->deadcrootc_transfer = max_deadcrootc * epc->woody_turnover;
			cs->deadcrootc = max_deadcrootc - cs->deadcrootc_transfer;
			ns->deadcrootn_transfer = cs->deadcrootc_transfer / epc->deadwood_cn;
			ns->deadcrootn = cs->deadcrootc / epc->deadwood_cn;

		}
	
		/* *****************************************************************************- */
		/* 4. use then penology array information to determine, for the first day of simulation, how many days of transfer and litterfall have already occurred for this year */
		
		if (phen->predays_transfer > 0)
		{
			prop_transfer = (double)phen->predays_transfer/(double)(phen->predays_transfer+phen->remdays_transfer);
			/* perform these transfers */
			transfer = prop_transfer * cs->leafc_transfer;
			cs->leafc          += transfer;
			cs->leafc_transfer -= transfer;
			transfer = prop_transfer * ns->leafn_transfer;
			ns->leafn          += transfer;
			ns->leafn_transfer -= transfer;
			transfer = prop_transfer * cs->frootc_transfer;
			cs->frootc          += transfer;
			cs->frootc_transfer -= transfer;
			transfer = prop_transfer * ns->frootn_transfer;
			ns->frootn          += transfer;
			ns->frootn_transfer -= transfer;
			/* yield simulation  */
			transfer = prop_transfer * cs->yieldc_transfer;
			cs->yieldc          += transfer;
			cs->yieldc_transfer -= transfer;
			transfer = prop_transfer * ns->yieldn_transfer;
			ns->yieldn          += transfer;
			ns->yieldn_transfer -= transfer;

			/* TREE-specific and NON-WOODY SPECIFIC fluxes */
			if (epc->woody)
			{
				transfer = prop_transfer * cs->livestemc_transfer;
				cs->livestemc          += transfer;
				cs->livestemc_transfer -= transfer;
				transfer = prop_transfer * ns->livestemn_transfer;
				ns->livestemn          += transfer;
				ns->livestemn_transfer -= transfer;
				transfer = prop_transfer * cs->deadstemc_transfer;
				cs->deadstemc          += transfer;
				cs->deadstemc_transfer -= transfer;
				transfer = prop_transfer * ns->deadstemn_transfer;
				ns->deadstemn          += transfer;
				ns->deadstemn_transfer -= transfer;
				transfer = prop_transfer * cs->livecrootc_transfer;
				cs->livecrootc          += transfer;
				cs->livecrootc_transfer -= transfer;
				transfer = prop_transfer * ns->livecrootn_transfer;
				ns->livecrootn          += transfer;
				ns->livecrootn_transfer -= transfer;
				transfer = prop_transfer * cs->deadcrootc_transfer;
				cs->deadcrootc          += transfer;
				cs->deadcrootc_transfer -= transfer;
				transfer = prop_transfer * ns->deadcrootn_transfer;
				ns->deadcrootn          += transfer;
				ns->deadcrootn_transfer -= transfer;
			}
			else
			{
				/* SOFT STEM SIMULATION of non-woody biomes */
				transfer = prop_transfer * cs->softstemc_transfer;
				cs->softstemc          += transfer;
				cs->softstemc_transfer -= transfer;
				transfer = prop_transfer * ns->softstemn_transfer;
				ns->softstemn          += transfer;
				ns->softstemn_transfer -= transfer;

			}
		
			/* only test for litterfall if there has already been some transfer growth this year */
			if (phen->predays_litfall > 0)
			{
				/* some litterfall has already occurred. in this case, just
				remove material from the displayed compartments, and don't
				bother with the transfer to litter compartments */
				prop_litfall = (double)phen->predays_litfall/(double)(phen->predays_litfall+phen->remdays_litfall);
				cs->leafc  -= prop_litfall * cs->leafc * epc->nonwoody_turnover;
				cs->frootc -= prop_litfall * cs->frootc * epc->nonwoody_turnover;
				cs->yieldc -= prop_litfall * cs->yieldc * epc->nonwoody_turnover;
				cs->softstemc -= prop_litfall * cs->softstemc * epc->nonwoody_turnover;
			}
		} /* end if transfer */
		/* add the growth respiration requirement for the first year's leaf and fine root growth from transfer pools to the gresp_transfer pool */
		cs->gresp_transfer = 0.0;
		cs->gresp_transfer += (cs->leafc_transfer + cs->frootc_transfer + cs->yieldc_transfer + cs->softstemc_transfer +
			                   cs->livestemc_transfer + cs->deadstemc_transfer + cs->livecrootc_transfer + cs->deadcrootc_transfer) * epc->GR_ratio; 
	

	}

	/*  if planting: transfer pools are deplenished (set to zero) on the first simulation day */ 
	if (PLT->PLT_num > 0 || ctrl->bareground_flag == 1)
	{
		cs->HRV_snkC += cs->leafc + cs->leafc_transfer;
		cs->HRV_snkC += cs->frootc + cs->frootc_transfer;
		cs->HRV_snkC += cs->yieldc + cs->yieldc_transfer;
		cs->HRV_snkC += cs->softstemc + cs->softstemc_transfer;
		cs->HRV_snkC += cs->gresp_transfer;
		cs->HRV_snkC += cs->STDBc_leaf;
		cs->HRV_snkC += cs->STDBc_froot;
		cs->HRV_snkC += cs->STDBc_yield;
		cs->HRV_snkC += cs->STDBc_softstem;

		cs->leafc = 0;
		cs->frootc = 0;
		cs->yieldc = 0;
		cs->softstemc = 0;
		cs->leafc_transfer = 0;
		cs->frootc_transfer = 0;
		cs->yieldc_transfer = 0;
		cs->softstemc_transfer = 0;
		cs->gresp_transfer = 0;
		cs->STDBc_leaf = 0;
		cs->STDBc_froot = 0;
		cs->STDBc_yield = 0;
		cs->STDBc_softstem = 0;

		ns->HRVsnk_N += ns->leafn + ns->leafn_transfer;
		ns->HRVsnk_N += ns->frootn + ns->frootn_transfer;
		ns->HRVsnk_N += ns->yieldn + ns->yieldn_transfer;
		ns->HRVsnk_N += ns->softstemn + ns->softstemn_transfer;
		ns->HRVsnk_N += ns->retransn;
		ns->HRVsnk_N += ns->STDBn_leaf;
		ns->HRVsnk_N += ns->STDBn_froot;
		ns->HRVsnk_N += ns->STDBn_yield;
		ns->HRVsnk_N += ns->STDBn_softstem;

		ns->leafn = 0;
		ns->frootn = 0;
		ns->yieldn = 0;
		ns->softstemn = 0;
		ns->leafn_transfer = 0;
		ns->frootn_transfer = 0;
		ns->yieldn_transfer = 0;
		ns->softstemn_transfer = 0;
		ns->retransn = 0;
		ns->STDBn_leaf = 0;
		ns->STDBn_froot = 0;
		ns->STDBn_yield = 0;
		ns->STDBn_softstem = 0;
	}
	
	/* set the initial rates of litterfall and live wood turnover */
	if (epc->evergreen)
	{
		/* leaf and fineroot litterfall rates */
		if (ctrl->spinup != 0)
		{
			epv->day_leafc_litfall_increment = cinit->max_leafc * epc->nonwoody_turnover / nDAYS_OF_YEAR;
			epv->day_frootc_litfall_increment = cinit->max_frootc * epc->nonwoody_turnover / nDAYS_OF_YEAR;
			epv->day_yield_litfall_increment = cinit->max_yieldc * epc->nonwoody_turnover / nDAYS_OF_YEAR;
			epv->day_softstemc_litfall_increment = cinit->max_softstemc * epc->nonwoody_turnover / nDAYS_OF_YEAR;
		}
		else
		{
			epv->day_leafc_litfall_increment = epv->annmax_leafc * epc->nonwoody_turnover / nDAYS_OF_YEAR;
			epv->day_frootc_litfall_increment = epv->annmax_frootc * epc->nonwoody_turnover / nDAYS_OF_YEAR;
			epv->day_yield_litfall_increment = epv->annmax_yieldc * epc->nonwoody_turnover / nDAYS_OF_YEAR;
			epv->day_softstemc_litfall_increment = epv->annmax_softstemc * epc->nonwoody_turnover / nDAYS_OF_YEAR;
		}
	}
	else
	{
		/* deciduous: reset the litterfall rates to 0.0 for the start of the
		next litterfall season */
		epv->day_leafc_litfall_increment = 0.0;
		epv->day_frootc_litfall_increment = 0.0;
		epv->day_yield_litfall_increment = 0.0;
		epv->day_softstemc_litfall_increment = 0.0;
	}

	epv->day_livestemc_turnover_increment = cs->livestemc * epc->woody_turnover / nDAYS_OF_YEAR;
	epv->day_livecrootc_turnover_increment = cs->livecrootc * epc->woody_turnover / nDAYS_OF_YEAR;


	/* in case of land-use change: agroecosystem to natural vegetation */
 	if (!ctrl->spinup && !PLT->PLT_num && ctrl->bareground_flag == 0)
	{
		if (!cs->leafc_transfer)
		{
			cs->leafc_transfer = 0.001;
			ns->leafn_transfer = cs->leafc_transfer/ epc->leaf_cn;
		}

		if (!cs->frootc_transfer)
		{
			cs->frootc_transfer = 0.001;
			ns->frootn_transfer = ns->frootn_transfer/ epc->leaf_cn;
		}
		
	}

	/* in case of land-use change: non-woody to woody and woody to non-woody */
	if (ctrl->spinup == 0)
	{
		if (epc->woody)
		{
			cs->softstemc = 0;
			cs->softstemc_storage = 0;
			cs->softstemc_transfer = 0;
			cs->HRV_snkC += cs->softstemc + cs->softstemc_storage + cs->softstemc_transfer;
			ns->softstemn = 0;
			ns->softstemn_storage = 0;
			ns->softstemn_transfer = 0;
			ns->HRVsnk_N += ns->softstemn + ns->softstemn_storage + ns->softstemn_transfer;
		}
		else
		{
			cs->livestemc = 0;
			cs->livestemc_storage = 0;
			cs->livestemc_transfer = 0;
			cs->livecrootc = 0;
			cs->livecrootc_storage = 0;
			cs->livecrootc_transfer = 0;
			cs->deadstemc = 0;
			cs->deadstemc_storage = 0;
			cs->deadstemc_transfer = 0;
			cs->deadcrootc = 0;
			cs->deadcrootc_storage = 0;
			cs->deadcrootc_transfer = 0;
			cs->HRV_snkC += cs->livestemc + cs->livestemc_storage + cs->livestemc_transfer +
				                  cs->livecrootc + cs->livecrootc_storage + cs->livecrootc_transfer +
								  cs->deadstemc + cs->deadstemc_storage + cs->deadstemc_transfer +
				                  cs->deadcrootc + cs->deadcrootc_storage + cs->deadcrootc_transfer;
			ns->livestemn = 0;
			ns->livestemn_storage = 0;
			ns->livestemn_transfer = 0;
			ns->livecrootn = 0;
			ns->livecrootn_storage = 0;
			ns->livecrootn_transfer = 0;
			ns->deadstemn = 0;
			ns->deadstemn_storage = 0;
			ns->deadstemn_transfer = 0;
			ns->deadcrootn = 0;
			ns->deadcrootn_storage = 0;
			ns->deadcrootn_transfer = 0;
			ns->HRVsnk_N += ns->livestemn + ns->livestemn_storage + ns->livestemn_transfer +
				                  ns->livecrootn + ns->livecrootn_storage + ns->livecrootn_transfer +
								  ns->deadstemn + ns->deadstemn_storage + ns->deadstemn_transfer +
				                  ns->deadcrootn + ns->deadcrootn_storage + ns->deadcrootn_transfer;
		}

	}

	/* initialization of dissolving coefficent array */
	soilInfo->dissolv_prop[0] = sprop->NH4_mobilen_prop;
	soilInfo->dissolv_prop[1] = NO3_mobilen_prop;
	soilInfo->dissolv_prop[2] = sprop->SOIL1dissolv_prop;
	soilInfo->dissolv_prop[3] = sprop->SOIL2dissolv_prop;
	soilInfo->dissolv_prop[4] = sprop->SOIL3dissolv_prop;
	soilInfo->dissolv_prop[5] = sprop->SOIL4dissolv_prop;
	soilInfo->dissolv_prop[6] = sprop->SOIL1dissolv_prop;
	soilInfo->dissolv_prop[7] = sprop->SOIL2dissolv_prop;
	soilInfo->dissolv_prop[8] = sprop->SOIL3dissolv_prop;
	soilInfo->dissolv_prop[9] = sprop->SOIL4dissolv_prop;

	
	

	/* call soil concentration calculation routine to calculate the concetration of soil (-1: all layers, NH4 -> content_soil*/
	if (!errorCode && check_soilcontent(-1, 0, sprop, cs, ns, soilInfo))
	{
		printf("ERROR in check_soilcontent.c for firstday.c\n");
		errorCode = 1;
	}

	/* hydroparams */
	if (!errorCode && multilayer_hydrolparams(sitec, sprop, ws, epv))
	{
		printf("\n");
		printf("ERROR in multilayer_hydrolparams.c for firstday.c\n");
		errorCode = 40701;
	}

	if (!errorCode && multilayer_rootDepth(epc, sprop, cs, sitec, epv))
	{
		printf("\n");
		printf("ERROR in multilayer_rootDepth.c from firstday.c\n");
		errorCode = 40702;
	}

	return (errorCode);
}
