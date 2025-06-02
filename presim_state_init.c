/* 
presim_state_init.c
Initialize water, carbon, and nitrogen state variables to 0.0 before
each simulation.

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGCMuSo v7.0.
Original code: Copyright 2000, Peter E. Thornton
Numerical Terradynamic Simulation Group, The University of Montana, USA
Modified code: Copyright 2022 D. Hidy [dori.hidy@gmail.com]
Hungarian Academy of Sciences, Hungary
See the website of Biome-BGCMuSo at http://nimbus.elte.hu/bbgc/ for documentation, model executable and example input files.
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
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

int presim_state_init(wstate_struct* ws, cstate_struct* cs, nstate_struct* ns, cinit_struct* cinit, control_struct* ctrl, soilInfo_struct* soilInfo)
{
	int errorCode=0;
	int layer, dm;
	
	/* initialization */
	ctrl->simyr = 0;
	ctrl->yday = 0;
	ctrl->plantyr = -1;
	ctrl->spinyears = 0;
	ctrl->month = 1;
	ctrl->day = 1;
	ctrl->limitEVP_flag = 0;
	ctrl->limitTRP_flag = 0;
	ctrl->limitMR_flag = 0;
	ctrl->limitDENIT_flag = 0;
	ctrl->limitSNSC_flag = 0;
	ctrl->limitleach_flag = 0;
	ctrl->limitdiffus_flag = 0;
	ctrl->pond_flag = 0;
	ctrl->noTRP_flag = 0;
	ctrl->grazingW_flag = 0;
	ctrl->condMOWerr_flag = 0;
	ctrl->condIRGerr_flag = 0;
	ctrl->prephen1_flag = 0;
	ctrl->prephen2_flag = 0;
	ctrl->bareground_flag = 0;
	ctrl->vegper_flag = 0;
	ctrl->south_shift = 0;
	ctrl->allocControl_flag = 0;
	ctrl->NaddSPINUP_flag = 0;

	ctrl->phenology_flag = 0;
	ctrl->transferGDD_flag = 0;
	ctrl->q10depend_flag = 0;
	ctrl->phtsyn_acclim_flag = 0;
	ctrl->resp_acclim_flag = 0;
	ctrl->CO2conduct_flag = 0;
	ctrl->STCM_flag = 0;
	ctrl->photosynt_flag = 0;
	ctrl->ET_flag = 0;
	ctrl->radiation_flag = 0;
	ctrl->soilstress_flag = 0;
	ctrl->interception_flag = 0;
	ctrl->MRdeficit_flag = 0;
	ctrl->Ksat_flag = 0;
	ctrl->firstsimday_flag = 1;

	for (layer = 0; layer < N_SOILLAYERS; layer++) ctrl->soiltype_array[layer] = 0;

	cinit->max_leafc = 0.0;
	cinit->max_frootc = 0.0;
	cinit->max_yieldc = 0.0;
	cinit->max_softstemc = 0.0;
	cinit->max_livestemc = 0.0;
	cinit->max_livecrootc = 0.0;
	
	ws->soilw_SUM = 0;
	ws->soilw_RZ = 0;
	ws->soilwAVAIL_RZ=0;
	ws->soilw_2m = 0;
	ws->pondw = 0;
	ws->snoww = 0;
	ws->canopyw = 0;
	ws->prcp_src = 0;
	ws->soilEVP_snk = 0;
	ws->snowSUBL_snk = 0;
	ws->EVPcanopyw_snk = 0;
	ws->TRP_snk = 0;
	ws->runoff_snk = 0;
	ws->pondEVP_snk = 0;
	ws->deeppercolation_snk = 0;
	ws->GWsrc_W = 0;
	ws->GWsnk_W = 0;
	ws->FLsrc_W = 0;
	ws->canopyw_THNsnk = 0;
	ws->canopyw_MOWsnk = 0;
	ws->canopyw_HRVsnk = 0;
	ws->canopyw_PLGsnk = 0;
	ws->canopyw_GRZsnk = 0;
	ws->IRGsrc_W = 0;
	ws->condIRGsrc = 0;
	ws->FRZsrc_W = 0;
	ws->GW_waterlogging = 0.0;
	ws->WbalanceERR = 0;
	ws->inW = 0;
	ws->outW = 0;
	ws->storeW = 0;
	ws->EVPsurface1cum = 0.0;
	ws->EVPsurface2cum = 0.0;
	cs->leafc = 0;
	cs->leafc_storage = 0;
	cs->leafc_transfer = 0;
	cs->frootc = 0;
	cs->frootc_storage = 0;
	cs->frootc_transfer = 0;
	cs->yieldc = 0;
	cs->yieldc_storage = 0;
	cs->yieldc_transfer = 0;
	cs->softstemc = 0;
	cs->softstemc_storage = 0;
	cs->softstemc_transfer = 0;
	cs->livestemc = 0;
	cs->livestemc_storage = 0;
	cs->livestemc_transfer = 0;
	cs->deadstemc = 0;
	cs->deadstemc_storage = 0;
	cs->deadstemc_transfer = 0;
	cs->livecrootc = 0;
	cs->livecrootc_storage = 0;
	cs->livecrootc_transfer = 0;
	cs->deadcrootc = 0;
	cs->deadcrootc_storage = 0;
	cs->deadcrootc_transfer = 0;
	cs->gresp_storage = 0;
	cs->gresp_transfer = 0;
	cs->litr1c_total = 0;
	cs->litr2c_total = 0;
	cs->litr3c_total = 0;
	cs->litr4c_total = 0;
	cs->cwdc_total = 0;
	cs->STDBc_leaf = 0;
	cs->STDBc_froot = 0;
	cs->STDBc_yield = 0;
	cs->STDBc_softstem = 0;
	cs->STDBc_above = 0;
	cs->STDBc_below = 0;
	cs->CTDBc_leaf = 0;
	cs->CTDBc_froot = 0;
	cs->CTDBc_yield = 0;
	cs->CTDBc_softstem = 0;
	cs->CTDBc_cstem = 0;
	cs->CTDBc_croot = 0;
	cs->CTDBc_above = 0;
	cs->CTDBc_below = 0;
	cs->soil1c_total = 0;
	cs->soil2c_total = 0;
	cs->soil3c_total = 0;
	cs->soil4c_total = 0;			
	cs->cpool = 0;
	cs->psnsun_src = 0;
	cs->psnshade_src = 0;
	cs->MRleaf_snk = 0;
	cs->GRleaf_snk = 0;
	cs->MRfroot_snk = 0;
	cs->GRfroot_snk = 0;
	cs->GRyield_snk = 0;
	cs->MRyield_snk = 0;
	cs->GRsoftstem_snk = 0;
	cs->MRsoftstem_snk = 0;
	cs->MRlivestem_snk = 0;
	cs->GRlivestem_snk = 0;
	cs->GRdeadstem_snk = 0;
	cs->MRlivecroot_snk = 0;
	cs->GRlivecroot_snk = 0;
	cs->GRdeadcroot_snk = 0;
	cs->HRlitr1_snk = 0;
	cs->HRlitr2_snk = 0;
	cs->HRlitr4_snk = 0;
	cs->HRsoil1_snk = 0;
	cs->HRsoil2_snk = 0;
	cs->HRsoil3_snk = 0;
	cs->HRsoil4_snk = 0;
	cs->Cdeepleach_snk = 0;
	cs->flowHSsnk_C = 0;
	cs->calc_flowHS = 0;
	cs->FIREsnk_C = 0;
	cs->SNSCsnk_C = 0;
	cs->PLTsrc_C = 0;
	cs->MULsrc_C = 0;
	cs->THN_transportC = 0;
	cs->HRV_snkC = 0;
	cs->MOW_snkC = 0;
	cs->GRZsnk_C = 0;
	cs->GRZsrc_C = 0;
	cs->FRZsrc_C = 0;
	cs->GWsrc_C = 0;
	cs->GWsnk_C = 0;
	cs->FLsrc_C = 0;
	cs->yieldC_HRV = 0.0;
	cs->frootC_HRV = 0.0;
	cs->vegCabove_HRV = 0.0;

	cs->CbalanceERR = 0;
	cs->CNratioERR = 0.0;
	cs->inC = 0;
	cs->outC = 0;
	cs->storeC = 0;

	ns->leafn = 0;
	ns->leafn_storage = 0;
	ns->leafn_transfer = 0;
	ns->frootn = 0;
	ns->frootn_storage = 0;
	ns->frootn_transfer = 0;
	ns->yieldn = 0;
	ns->yieldn_storage = 0;
	ns->yieldn_transfer = 0;
	ns->softstemn = 0;
	ns->softstemn_storage = 0;
	ns->softstemn_transfer = 0;
	ns->livestemn = 0;
	ns->livestemn_storage = 0;
	ns->livestemn_transfer = 0;
	ns->deadstemn = 0;
	ns->deadstemn_storage = 0;
	ns->deadstemn_transfer = 0;
	ns->livecrootn = 0;
	ns->livecrootn_storage = 0;
	ns->livecrootn_transfer = 0;
	ns->deadcrootn = 0;
	ns->deadcrootn_storage = 0;
	ns->deadcrootn_transfer = 0;
	ns->npool = 0;
	ns->litr1n_total = 0;
	ns->litr2n_total = 0;
	ns->litr3n_total = 0;
	ns->litr4n_total = 0;
	ns->cwdn_total = 0;
	ns->STDBn_leaf = 0;
	ns->STDBn_froot = 0;
	ns->STDBn_yield = 0;
	ns->STDBn_softstem = 0;
	ns->STDBn_above = 0;
	ns->STDBn_below = 0;
	ns->CTDBn_leaf = 0;
	ns->CTDBn_froot = 0;
	ns->CTDBn_yield = 0;
	ns->CTDBn_softstem = 0;
	ns->CTDBn_cstem = 0;
	ns->CTDBn_croot = 0;
	ns->CTDBn_above = 0;
	ns->CTDBn_below = 0;
	ns->soil1n_total = 0;
	ns->soil2n_total = 0;
	ns->soil3n_total = 0;
	ns->soil4n_total = 0;
	ns->retransn = 0;
	ns->NH4_total = 0;
	ns->NO3_total = 0;
	ns->Nfix_src = 0;
	ns->Ndep_src = 0;
	ns->Ndeepleach_snk = 0;
	ns->Nvol_snk = 0;
	ns->Nprec_snk = 0;
	ns->FIREsnk_N = 0;
	ns->SNSCsnk_N = 0;
	ns->FRZsrc_N = 0;
	ns->PLTsrc_N = 0;
	ns->MULsrc_N = 0;
	ns->CWEsnk_N = 0;
	ns->FLsrc_N = 0;
	ns->GWsrc_N = 0;
	ns->GWsnk_N = 0;
	ns->THNsnk_N = 0;
	ns->MOWsnk_N = 0;
	ns->HRVsnk_N = 0;
	ns->GRZsnk_N = 0;
	ns->GRZsrc_N = 0;
	
	ns->SPINUPsrc = 0;
	ns->NbalanceERR = 0;
	ns->inN = 0;
	ns->outN = 0;
	ns->storeN = 0;


	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		ws->soilw[layer] = 0;
		ws->soilw_pre[layer] = 0;
		ws->soilwAVAIL[layer] = 0;
		cs->cwdc[layer] = 0;
		cs->litr1c[layer] = 0;
		cs->litr2c[layer] = 0;
		cs->litr3c[layer] = 0;
		cs->litr4c[layer] = 0;
		cs->litrC[layer] = 0;
		cs->soilC[layer] = 0;
		ns->cwdn[layer] = 0;
		ns->litr1n[layer] = 0;
		ns->litr2n[layer] = 0;
		ns->litr3n[layer] = 0;
		ns->litr4n[layer] = 0;
		ns->litrN[layer] = 0;
		ns->soil1n[layer] = 0;
		ns->soil2n[layer] = 0;
		ns->soil3n[layer] = 0;
		ns->soil4n[layer] = 0;
		ns->soilN[layer] = 0;
		ns->NH4[layer] = 0;
		ns->NO3[layer] = 0;	
	}

	for (dm = 0; dm < N_DISSOLVMATER; dm++)
	{
		soilInfo->dissolv_prop[dm] = 0;
		soilInfo->GWconc[dm] = 0;
		soilInfo->FLconc[dm] = 0;
		soilInfo->content_NORMgw[dm] = 0;
		soilInfo->content_CAPILgw[dm] = 0;
		soilInfo->content_SATgw[dm] = 0;
		soilInfo->dismatLeach_NORMgw[dm] = 0;
		soilInfo->dismatLeach_NORMcf[dm] = 0;
		soilInfo->dismatGWmovchange[dm] = 0;
		soilInfo->dismatGWecofunc_NORM[dm] = 0;
		soilInfo->dismatGWecofunc_CAPIL[dm] = 0;
		soilInfo->dismatGWdecomp_CAPIL[dm] = 0;
		soilInfo->dismatGWdecomp_NORM[dm] = 0;


		for (layer = 0; layer < N_SOILLAYERS; layer++)
		{
			soilInfo->content_soil[dm][layer] = 0;
			soilInfo->dismatLeach[dm][layer] = 0;
			soilInfo->dismatLeach_percolDiffus[dm][layer] = 0;
			soilInfo->dismatGWdischarge[dm][layer] = 0;
			soilInfo->dismatGWrecharge[dm][layer] = 0;
			soilInfo->dismatGWecofunc[dm][layer] = 0;
			soilInfo->dismatGWdecomp[dm][layer] = 0;
			soilInfo->dismatGWfertil[dm][layer] = 0;
			soilInfo->dismatUNSATecofunc[dm][layer] = 0;
			soilInfo->dismatUNSATdecomp[dm][layer] = 0;
			soilInfo->dismatUNSATfertil[dm][layer] = 0;
			soilInfo->dismatTOTALecofunc[dm][layer] = 0;
			soilInfo->dismatTOTALdecomp[dm][layer] = 0;
			soilInfo->dismatTOTALfertil[dm][layer] = 0;

			if (dm < N_DISSOLVorgN)
			{
				soilInfo->FRZ_to_litrN[dm][layer] = 0;
				soilInfo->FRZ_to_litrC[dm][layer] = 0;
			}

		}
	}


	return(errorCode);
}
	
