/* 
summary.c
summary variables for potential output

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
#include <string.h>
#include <math.h>
#include <malloc.h>
#include "ini.h"
#include "bgc_struct.h"
#include "bgc_func.h"
#include "bgc_constants.h"

int cnw_summary(const epconst_struct* epc, const siteconst_struct* sitec, const soilprop_struct* sprop, const metvar_struct* metv,
	const cstate_struct* cs, const cflux_struct* cf, const nstate_struct* ns, const nflux_struct* nf, const wflux_struct* wf, const soilInfo_struct* soilInfo,
	epvar_struct* epv, summary_struct* summary)
{
	int errorCode = 0;
	int layer, GWlayer, dm, N_SOILLAYERS30cm, N_SOILLAYERS10cm;
	double GPP, MR, GR, HR, TR, AR, fire;
	double SR; /* calculating soil respiration */
	double NPP, NEP, NEE, NBP, disturb_loss, disturb_gain, BD_top30, BD_top10, BD_act, SOC_top10, prop_to_percent;
	double NPPabove_w, NPPbelow_w, NPPabove_nw, NPPbelow_nw;
	double Closs_THN_w, Closs_THN_nw, Closs_MOW, Closs_HRV, yieldC_HRV, Closs_PLG, Closs_PLT, Closs_GRZ, Cplus_PLT, Cplus_FRZ, Cplus_GRZ;
	double Nplus_PLT, Nloss_PLT, Nplus_GRZ, Nplus_FRZ;
	double Closs_SNSC, STDB_to_litr, CTDB_to_litr;
	double GRabove_nw, GRbelow_nw, GRabove_w, GRbelow_w, MRabove_nw, MRbelow_nw, MRabove_w, MRbelow_w, MRdef_fluxes, diff;
	double vegN_above, vegN_below, unSATprop;


	summary->leafCN = summary->frootCN = summary->yieldN = summary->softstemCN = 0;
	NPP = NPPabove_w = NPPbelow_w = NPPabove_nw = NPPbelow_nw = 0;

	/* actual phenological phase */
	int ap = (int)epv->n_actphen - 1;

	GWlayer = (int)sprop->GWlayer;

	N_SOILLAYERS30cm = 3;
	N_SOILLAYERS10cm = 2;

	/************************************************************************************************************************************/
	/* 1. summarize meteorological and water variables */


	summary->annprcp += metv->prcp;
	summary->anntavg += metv->Tavg / nDAYS_OF_YEAR;


	/************************************************************************************************************************************/
	/* 2. rootzone leaching variables: kg to g */

	summary->sminNleach_RZmax = 0;
	summary->DONleach_RZmax = 0;
	summary->DOCleach_RZmax = 0;
	for (dm = 0; dm < N_DISSOLVMATER; dm++)
	{
		if (dm < N_DISSOLVinorgN)
			summary->sminNleach_RZmax += soilInfo->content_soil[dm][epv->n_maxrootlayers - 1];
		else
		{
			if (dm < N_DISSOLVN)
				summary->DONleach_RZmax += soilInfo->content_soil[dm][epv->n_maxrootlayers - 1];
			else
				summary->DOCleach_RZmax += soilInfo->content_soil[dm][epv->n_maxrootlayers - 1];
		}

	}


	summary->Wleach_RZmax = wf->soilwFlux[epv->n_maxrootlayers - 1];


	summary->cumRunoff += wf->prcp_to_runoff + wf->pondw_to_runoff;
	summary->cumWleach_RZmax += summary->Wleach_RZmax;
	summary->cumNleach_RZmax += summary->sminNleach_RZmax;

	/************************************************************************************************************************************/
	/* 3. summarize carbon and nitrogen stocks */


	/*  biomass C (live+dead) without NSC */
	summary->LDaboveC_nw = cs->leafc + cs->yieldc + cs->softstemc +
		cs->STDBc_leaf + cs->STDBc_yield + cs->STDBc_softstem;
	summary->LDaboveC_w = cs->livestemc + cs->deadstemc;

	summary->LDbelowC_nw = cs->frootc + cs->STDBc_froot;
	summary->LDbelowC_w = cs->livecrootc + cs->deadcrootc;

	/* biomass C (live+dead) with NSC */
	summary->LDaboveCwithNSC_nw = cs->leafc + cs->yieldc + cs->softstemc +
		cs->STDBc_above +
		cs->leafc_transfer + cs->leafc_storage +
		cs->yieldc_storage + cs->yieldc_transfer + cs->softstemc_storage + cs->softstemc_transfer +
		cs->gresp_storage + cs->gresp_transfer;
	summary->LDaboveCwithNSC_w = cs->livestemc + cs->deadstemc +
		cs->livestemc_storage + cs->livestemc_transfer +
		cs->deadstemc_storage + cs->deadstemc_transfer;

	summary->LDbelowCwithNSC_nw = cs->frootc +
		cs->STDBc_below +
		cs->frootc_storage + cs->frootc_transfer;

	summary->LDbelowCwithNSC_w = cs->livecrootc + cs->deadcrootc +
		cs->livecrootc_storage + cs->livecrootc_transfer +
		cs->deadcrootc_storage + cs->deadcrootc_transfer;


	/* living biomass C  */
	summary->LaboveC_nw = cs->leafc + cs->yieldc + cs->softstemc;
	summary->LaboveC_w = cs->livestemc;

	summary->LbelowC_nw = cs->frootc;
	summary->LbelowC_w = cs->deadcrootc;

	/* living biomass C with NSC */
	summary->LaboveCwithNSC_nw = cs->leafc + cs->yieldc + cs->softstemc +
		cs->leafc_transfer + cs->leafc_storage +
		cs->yieldc_storage + cs->yieldc_transfer + cs->softstemc_storage + cs->softstemc_transfer +
		cs->gresp_storage + cs->gresp_transfer;
	summary->LaboveCwithNSC_w = cs->livestemc + cs->livestemc_storage + cs->livestemc_transfer;

	summary->LbelowCwithNSC_nw = cs->frootc + cs->frootc_storage + cs->frootc_transfer;
	summary->LbelowCwithNSC_w = cs->deadstemc + cs->deadstemc_storage + cs->deadstemc_transfer;


	/* dead biomass C */
	summary->DaboveC_nw = cs->STDBc_leaf + cs->STDBc_yield + cs->STDBc_softstem;
	summary->DaboveC_w = cs->deadstemc;

	summary->DbelowC_nw = cs->STDBc_froot;
	summary->DbelowC_w = cs->deadcrootc;

	/* dead biomass C with NSC  */
	summary->DaboveCwithNSC_nw = cs->STDBc_leaf + cs->STDBc_yield + cs->STDBc_softstem;
	summary->DaboveCwithNSC_w = cs->deadstemc + cs->deadstemc_storage + cs->deadstemc_transfer;

	summary->DbelowCwithNSC_nw = cs->STDBc_froot;
	summary->DbelowCwithNSC_w = cs->deadstemc + cs->deadstemc_storage + cs->deadstemc_transfer;

	/* living SC and NSC */
	summary->livingSC = cs->leafc + cs->yieldc + cs->softstemc + cs->frootc + cs->livestemc + cs->deadstemc;
	summary->livingNSC = cs->leafc_storage + cs->leafc_transfer +
		cs->gresp_storage + cs->gresp_transfer +
		cs->yieldc_storage + cs->yieldc_transfer +
		cs->softstemc_storage + cs->softstemc_transfer +
		cs->livestemc_storage + cs->livestemc_transfer +
		cs->deadstemc_storage + cs->deadstemc_transfer +
		cs->frootc_storage + cs->frootc_transfer +
		cs->livecrootc_storage + cs->livecrootc_transfer +
		cs->deadcrootc_storage + cs->deadcrootc_transfer;;

	/* living and total, above- and belowground biomass (C+N) */
	summary->livingBIOMabove = cs->leafc + cs->leafc_storage + cs->leafc_transfer +
		cs->gresp_storage + cs->gresp_transfer +
		cs->yieldc + cs->yieldc_storage + cs->yieldc_transfer +
		cs->softstemc + cs->softstemc_storage + cs->softstemc_transfer +
		cs->livestemc + cs->livestemc_storage + cs->livestemc_transfer +
		cs->deadstemc + cs->deadstemc_storage + cs->deadstemc_transfer +
		ns->leafn + ns->leafn_storage + ns->leafn_transfer +
		ns->retransn +
		ns->yieldn + ns->yieldn_storage + ns->yieldn_transfer +
		ns->softstemn + ns->softstemn_storage + ns->softstemn_transfer +
		ns->livestemn + ns->livestemn_storage + ns->livestemn_transfer +
		ns->deadstemn + ns->deadstemn_storage + ns->deadstemn_transfer;
	summary->livingBIOMbelow = cs->frootc + cs->frootc_storage + cs->frootc_transfer +
		cs->livecrootc + cs->livecrootc_storage + cs->livecrootc_transfer +
		cs->deadcrootc + cs->deadcrootc_storage + cs->deadcrootc_transfer +
		ns->frootn + ns->frootn_storage + ns->frootn_transfer +
		ns->livecrootn + ns->livecrootn_storage + ns->livecrootn_transfer +
		ns->deadcrootn + ns->deadcrootn_storage + ns->deadcrootn_transfer;
	summary->BIOMaboveSUM = summary->livingBIOMabove + cs->STDBc_above + ns->STDBn_above;

	summary->BIOMbelowSUM = summary->livingBIOMbelow + cs->STDBc_below + ns->STDBn_below;


	if (summary->livingBIOMabove > summary->annmax_livingBIOMabove) summary->annmax_livingBIOMabove = summary->livingBIOMabove;
	if (summary->livingBIOMbelow > summary->annmax_livingBIOMbelow) summary->annmax_livingBIOMbelow = summary->livingBIOMbelow;
	if (summary->BIOMaboveSUM > summary->annmax_BIOMaboveSUM)  summary->annmax_BIOMaboveSUM = summary->BIOMaboveSUM;
	if (summary->BIOMbelowSUM > summary->annmax_BIOMbelowSUM)  summary->annmax_BIOMbelowSUM = summary->BIOMbelowSUM;

	summary->vegC = cs->leafc + cs->leafc_storage + cs->leafc_transfer +
		cs->frootc + cs->frootc_storage + cs->frootc_transfer +
		cs->yieldc + cs->yieldc_storage + cs->yieldc_transfer +
		cs->softstemc + cs->softstemc_storage + cs->softstemc_transfer +
		cs->livestemc + cs->livestemc_storage + cs->livestemc_transfer +
		cs->deadstemc + cs->deadstemc_storage + cs->deadstemc_transfer +
		cs->livecrootc + cs->livecrootc_storage + cs->livecrootc_transfer +
		cs->deadcrootc + cs->deadcrootc_storage + cs->deadcrootc_transfer +
		cs->gresp_storage + cs->gresp_transfer +
		cs->STDBc_above + cs->STDBc_below +
		cs->CTDBc_above + cs->CTDBc_below;

	summary->vegN = ns->leafn + ns->leafn_storage + ns->leafn_transfer +
		ns->frootn + ns->frootn_storage + ns->frootn_transfer +
		ns->yieldn + ns->yieldn_storage + ns->yieldn_transfer +
		ns->softstemn + ns->softstemn_storage + ns->softstemn_transfer +
		ns->livestemn + ns->livestemn_storage + ns->livestemn_transfer +
		ns->deadstemn + ns->deadstemn_storage + ns->deadstemn_transfer +
		ns->livecrootn + ns->livecrootn_storage + ns->livecrootn_transfer +
		ns->deadcrootn + ns->deadcrootn_storage + ns->deadcrootn_transfer +
		ns->retransn +
		ns->STDBn_above + ns->STDBn_below +
		ns->CTDBn_above + ns->CTDBn_below;

	summary->leafc_LandD = cs->leafc + cs->STDBc_leaf;
	summary->frootc_LandD = cs->frootc + cs->STDBc_froot;
	summary->yield_LandD = cs->yieldc + cs->STDBc_yield;
	summary->softstemc_LandD = cs->softstemc + cs->STDBc_softstem;

	summary->leafDM = summary->leafc_LandD / epc->leafC_DM;
	summary->frootDM = summary->frootc_LandD / epc->frootC_DM;
	summary->yieldDM = summary->yield_LandD / epc->yield_DM;
	summary->softstemDM = summary->softstemc_LandD / epc->softstemC_DM;


	if (ns->leafn + ns->STDBn_leaf != 0)         summary->leafCN = (cs->leafc + cs->STDBc_leaf) / (ns->leafn + ns->STDBn_leaf);
	if (ns->frootn + ns->STDBn_froot != 0)       summary->frootCN = (cs->frootc + cs->STDBc_froot) / (ns->frootn + ns->STDBn_froot);
	if (ns->yieldn + ns->STDBn_yield != 0)       summary->yieldN = (cs->yieldc + cs->STDBc_yield) / (ns->yieldn + ns->STDBn_yield);
	if (ns->softstemn + ns->STDBn_softstem != 0) summary->softstemCN = (cs->softstemc + cs->STDBc_softstem) / (ns->softstemn + ns->STDBn_softstem);

	summary->yieldDM_HRV = cs->yieldC_HRV / epc->yield_DM;

	summary->leaflitrDM = (cs->litr1c_total + cs->litr2c_total + cs->litr3c_total + cs->litr4c_total) / epc->leaflitrC_DM;
	summary->livewoodDM = (cs->livestemc + cs->livecrootc) / epc->livewoodC_DM;
	summary->deadwoodDM = (cs->deadstemc + cs->deadcrootc) / epc->deadwoodC_DM;



	summary->litrC_total = cs->litr1c_total + cs->litr2c_total + cs->litr3c_total + cs->litr4c_total;
	summary->litrN_total = ns->litr1n_total + ns->litr2n_total + ns->litr3n_total + ns->litr4n_total;

	summary->litrCwdC_total = cs->litr1c_total + cs->litr2c_total + cs->litr3c_total + cs->litr4c_total + cs->cwdc_total;
	summary->litrCwdN_total = ns->litr1n_total + ns->litr2n_total + ns->litr3n_total + ns->litr4n_total + ns->cwdn_total;;

	summary->soilC_total = cs->soil1c_total + cs->soil2c_total + cs->soil3c_total + cs->soil4c_total;
	summary->soilN_total = ns->soil1n_total + ns->soil2n_total + ns->soil3n_total + ns->soil4n_total;

	if (summary->litrN_total)
		summary->CNlitr_total = summary->litrC_total / summary->litrN_total;
	else
		summary->CNlitr_total = 0;

	if (summary->soilN_total)
		summary->CNsoil_total = summary->soilC_total / summary->soilN_total;
	else
		summary->CNsoil_total = 0;

	summary->sminN_total = ns->NH4_total + ns->NO3_total;
	summary->sminNdissolv_total = ns->NH4_total * sprop->NH4_mobilen_prop + ns->NO3_total * NO3_mobilen_prop;

	vegN_above = ns->leafn + ns->leafn_storage + ns->leafn_transfer +
		ns->yieldn + ns->yieldn_storage + ns->yieldn_transfer +
		ns->softstemn + ns->softstemn_storage + ns->softstemn_transfer +
		ns->livestemn + ns->livestemn_storage + ns->livestemn_transfer +
		ns->deadstemn + ns->deadstemn_storage + ns->deadstemn_transfer +
		ns->retransn +
		ns->STDBn_above + ns->CTDBn_above +
		ns->npool;

	vegN_below = ns->frootn + ns->frootn_storage + ns->frootn_transfer +
		ns->livecrootn + ns->livecrootn_storage + ns->livecrootn_transfer +
		ns->deadcrootn + ns->deadcrootn_storage + ns->deadcrootn_transfer +
		ns->STDBn_below + ns->CTDBn_below;


	summary->orgN = vegN_above + vegN_below + summary->litrN_total + summary->soilN_total;


	summary->cumMRdeficit_NSC += cf->NSCnw_to_MRdef + cf->NSCw_to_MRdef;
	summary->cumMRdeficit_SC += cf->SCnw_to_MRdef + cf->SCw_to_MRdef;

	summary->totalC = summary->vegC + summary->litrC_total + summary->soilC_total + cs->cwdc_total;
	summary->totalN = summary->vegN + summary->litrN_total + summary->soilN_total + ns->cwdn_total + summary->sminN_total;

	/************************************************************************************************************************************/
	/* carbon and nitrogen content of top soil layer (10 cm layer depth):
	   kg (C or N)/m2 -> g (C or N) / kg soil: kgC/m2 = kgCN/0.1m3 = 10 * kgCN/m3 */


	prop_to_percent = 100;

	summary->NH4_RZmax = 0;
	summary->NO3_RZmax = 0;
	summary->sminN_RZmax = 0;
	summary->soilC_RZmax = 0;
	summary->soilN_RZmax = 0;
	summary->litrC_RZmax = 0;
	summary->litrN_RZmax = 0;
	summary->sminNdissolv_RZmax = 0;

	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		/* NH4: kgN/m2; BD: g/cm3 -> kg/m3: *10-3; ppm: *1000000 */
		summary->NH4_ppm[layer] = (ns->NH4[layer] / (sprop->BD[layer] * g_per_cm3_to_kg_per_m3 * sitec->soillayer_thickness[layer])) * multi_ppm;
		summary->NO3_ppm[layer] = (ns->NO3[layer] / (sprop->BD[layer] * g_per_cm3_to_kg_per_m3 * sitec->soillayer_thickness[layer])) * multi_ppm;

		summary->sminNdissolv[layer] = ns->NH4[layer] * sprop->NH4_mobilen_prop + ns->NO3[layer] * NO3_mobilen_prop;

		if (layer < epv->n_maxrootlayers)
		{
			summary->sminN_RZmax += (ns->NH4[layer] + ns->NO3[layer]);
			summary->NO3_RZmax += ns->NO3[layer];
			summary->NH4_RZmax += ns->NH4[layer];
			summary->sminNdissolv_RZmax += (ns->NH4[layer] * sprop->NH4_mobilen_prop + ns->NO3[layer]);
			summary->soilC_RZmax += (cs->soilC[layer]);
			summary->soilN_RZmax += (ns->soilN[layer]);
			summary->litrC_RZmax += (cs->litrC[layer]);
			summary->litrN_RZmax += (ns->litrN[layer]);
		}
		BD_act = sprop->BD[layer] * g_per_cm3_to_kg_per_m3 * sitec->soillayer_thickness[layer];
		summary->SOCpercent[layer] = (cs->soil1c[layer] + cs->soil2c[layer] + cs->soil3c[layer] + cs->soil4c[layer]) / BD_act * prop_to_percent;
	}

	/* g/cm3 to kg/m2 */
	BD_top30 = 0;
	BD_top10 = 0;
	SOC_top10 = 0;
	for (layer = 0; layer < N_SOILLAYERS30cm; layer++)
	{
		if (layer < N_SOILLAYERS10cm)
		{ 
			BD_top30 += sprop->BD[layer] * sitec->soillayer_thickness[layer] * g_per_cm3_to_kg_per_m3;
			summary->SOC_top30 += cs->soil1c[layer] + cs->soil2c[layer] + cs->soil3c[layer] + cs->soil4c[layer];
			summary->SON_top30 += ns->soil1n[layer] + ns->soil2n[layer] + ns->soil3n[layer] + ns->soil4n[layer];
			summary->NH4dissolv_top30ppm += (summary->NH4_ppm[layer] * sitec->soillayer_thickness[layer] / sitec->soillayer_depth[N_SOILLAYERS30cm-1]) * sprop->NH4_mobilen_prop;
			summary->NO3dissolv_top30ppm += summary->NO3_ppm[layer] * sitec->soillayer_thickness[layer] / sitec->soillayer_depth[N_SOILLAYERS30cm-1];

			BD_top10 += sprop->BD[layer] * sitec->soillayer_thickness[layer] * g_per_cm3_to_kg_per_m3;
			SOC_top10 += cs->soil1c[layer] + cs->soil2c[layer] + cs->soil3c[layer] + cs->soil4c[layer];
		}
		else
		{
			BD_top30 += sprop->BD[layer] * sitec->soillayer_thickness[layer] * g_per_cm3_to_kg_per_m3;
			summary->SOC_top30 += cs->soil1c[layer] + cs->soil2c[layer] + cs->soil3c[layer] + cs->soil4c[layer];
			summary->SON_top30 += ns->soil1n[layer] + ns->soil2n[layer] + ns->soil3n[layer] + ns->soil4n[layer];
			summary->NH4dissolv_top30ppm += (summary->NH4_ppm[layer] * sitec->soillayer_thickness[layer] / sitec->soillayer_depth[N_SOILLAYERS30cm - 1]) * sprop->NH4_mobilen_prop;
			summary->NO3dissolv_top30ppm += summary->NO3_ppm[layer] * sitec->soillayer_thickness[layer] / sitec->soillayer_depth[N_SOILLAYERS30cm - 1];
		}
	}

	summary->SOCpercent_top10 = SOC_top10 / BD_top10 * prop_to_percent;
	summary->SOCpercent_top30 = summary->SOC_top30 / BD_top30 * prop_to_percent;
	summary->SONpercent_top30 = summary->SON_top30 / BD_top30 * prop_to_percent;
	summary->sminN_top30avail = summary->NO3dissolv_top30ppm + summary->NH4dissolv_top30ppm;

	summary->N2Oflux_total = nf->N2OfluxNITRIF_total + nf->N2OfluxDENITR_total + nf->N2OfluxGRZ + nf->N2OfluxFRZ;
	summary->CH4flux_total = cf->CH4flux_animal + cf->CH4flux_manure + cf->CH4flux_soil;

	/* GWP: CO2-eq from N2O and CH4 mass */
	summary->N2OfluxCeq = (summary->N2Oflux_total * N_to_N2O * GWP_N2O) * CO2_to_C;
	summary->CH4fluxCeq = (summary->CH4flux_total * C_to_CH4 * GWP_CH4) * CO2_to_C;

	/************************************************************************************************************************************/
	/* 5. calculate daily fluxes (GPP, NPP, NEP, MR, GR, HR) positive for net growth: NPP = gross PSN - Maintenance Resp - growth Resp */

	GPP = cf->psnsun_to_cpool + cf->psnshade_to_cpool;

	MRabove_nw = cf->leaf_day_MR + cf->leaf_night_MR + cf->yield_MR + cf->softstem_MR;
	MRbelow_nw = cf->froot_MR;
	MRabove_w = cf->livestem_MR;
	MRbelow_w = cf->livecroot_MR;

	MR = MRabove_nw + MRbelow_nw + MRabove_w + MRbelow_w;


	GRabove_nw = cf->cpool_leaf_GR + cf->cpool_leaf_storage_GR + cf->transfer_leaf_GR +
		cf->cpool_yield_GR + cf->cpool_yieldc_storage_GR + cf->transfer_yield_GR +
		cf->cpool_softstem_GR + cf->cpool_softstem_storage_GR + cf->transfer_softstem_GR;
	GRbelow_nw = cf->cpool_froot_GR + cf->cpool_froot_storage_GR + cf->transfer_froot_GR;
	GRabove_w = cf->cpool_livestem_GR + cf->cpool_livestem_storage_GR + cf->transfer_livestem_GR +
		cf->cpool_deadstem_GR + cf->cpool_deadstem_storage_GR + cf->transfer_deadstem_GR;
	GRbelow_w = cf->cpool_livecroot_GR + cf->cpool_livecroot_storage_GR + cf->transfer_livecroot_GR +
		cf->cpool_deadcroot_GR + cf->cpool_deadcroot_storage_GR + cf->transfer_deadcroot_GR;

	GR = GRabove_nw + GRbelow_nw + GRabove_w + GRbelow_w;


	NPP = GPP - MR - GR;

	/* NPP-partitioning */

	if (NPP > CRIT_PREC)
	{
		MRdef_fluxes = cf->NSCnw_to_MRdef + cf->SCnw_to_MRdef + cf->NSCw_to_MRdef + cf->SCw_to_MRdef + cf->cpool_to_MRdef;

		NPPabove_w = GPP * (epc->alloc_deadstemc[ap] + epc->alloc_livestemc[ap])
			- cf->livestem_MR - GRabove_w;

		NPPbelow_w = GPP * (epc->alloc_deadcrootc[ap] + epc->alloc_livecrootc[ap])
			- cf->livecroot_MR - GRbelow_w;

		NPPabove_nw = GPP * (epc->alloc_leafc[ap] + epc->alloc_softstemc[ap] + epc->alloc_yield[ap])
			- (cf->leaf_day_MR + cf->leaf_night_MR + cf->yield_MR + cf->softstem_MR) - GRabove_nw;

		NPPbelow_nw = GPP * (epc->alloc_frootc[ap]) - cf->froot_MR - GRbelow_nw;

		MR = MRabove_nw + MRbelow_nw + MRabove_w + MRbelow_w;

		diff = NPP - NPPabove_w - NPPbelow_w - NPPabove_nw - NPPbelow_nw;

		/* MRdef-fluxes can affect NPP partitioning calculation */
		if (fabs(diff) > CRIT_PREC)
		{
			if (diff > 0 && fabs(diff - MRdef_fluxes) < CRIT_PREC)
			{
				NPPabove_w += diff * (MRabove_w / MR);
				NPPbelow_w += diff * (MRbelow_w / MR);
				NPPabove_nw += diff * (MRabove_nw / MR);
				NPPbelow_nw += diff * (MRbelow_nw / MR);
			}
			else
			{
				printf("ERROR in NPP calculation (summary.c)\n");
				errorCode = 1;
			}
		}
	}
	else
		NPP = 0;


	/* heterotrophic respiration  */
	HR = 0;
	summary->litr1HR_total = 0;
	summary->litr2HR_total = 0;
	summary->litr4HR_total = 0;
	summary->soil1HR_total = 0;
	summary->soil2HR_total = 0;
	summary->soil3HR_total = 0;
	summary->soil4HR_total = 0;

	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		summary->litr1HR_total += cf->litr1_hr[layer];
		summary->litr2HR_total += cf->litr2_hr[layer];
		summary->litr4HR_total += cf->litr4_hr[layer];
		summary->soil1HR_total += cf->soil1_hr[layer];
		summary->soil2HR_total += cf->soil2_hr[layer];
		summary->soil3HR_total += cf->soil3_hr[layer];
		summary->soil4HR_total += cf->soil4_hr[layer];

		HR += cf->litr1_hr[layer] + cf->litr2_hr[layer] + cf->litr4_hr[layer] +
			cf->soil1_hr[layer] + cf->soil2_hr[layer] + cf->soil3_hr[layer] + cf->soil4_hr[layer];
	}

	TR = MR + GR + HR;

	AR = MR + GR;

	NEP = NPP - HR;

	/* soil respiration */
	SR = cf->froot_MR + cf->livecroot_MR +
		cf->cpool_froot_GR + cf->cpool_froot_storage_GR + cf->transfer_froot_GR +
		cf->cpool_livecroot_GR + cf->cpool_livecroot_storage_GR + cf->transfer_livecroot_GR +
		cf->cpool_deadcroot_GR + cf->cpool_deadcroot_storage_GR + cf->transfer_deadcroot_GR +
		HR;


	/* calculate daily NEE, positive for net sink: NEE = NEP - fire losses */
	fire = cf->m_leafc_to_fire + cf->m_leafc_storage_to_fire + cf->m_leafc_transfer_to_fire +
		cf->m_frootc_to_fire + cf->m_frootc_storage_to_fire + cf->m_frootc_transfer_to_fire +
		cf->m_livestemc_to_fire + cf->m_livestemc_storage_to_fire + cf->m_livestemc_transfer_to_fire +
		cf->m_deadstemc_to_fire + cf->m_deadstemc_storage_to_fire + cf->m_deadstemc_transfer_to_fire +
		cf->m_livecrootc_to_fire + cf->m_livecrootc_storage_to_fire + cf->m_livecrootc_transfer_to_fire +
		cf->m_deadcrootc_to_fire + cf->m_deadcrootc_storage_to_fire + cf->m_deadcrootc_transfer_to_fire +
		cf->m_gresp_storage_to_fire + cf->m_gresp_transfer_to_fire +
		cf->m_litr1c_to_fire_total + cf->m_litr2c_to_fire_total + cf->m_litr3c_to_fire_total + cf->m_litr4c_to_fire_total +
		cf->m_cwdc_to_fire_total +
		/* yield simulation */
		cf->m_yieldc_to_fire + cf->m_yieldc_storage_to_fire + cf->m_yieldc_transfer_to_fire +
		/* softstem simulation */
		cf->m_softstemc_to_fire + cf->m_softstemc_storage_to_fire + cf->m_softstemc_transfer_to_fire;
	/* NEE is positive if ecosystem is net source and negative if it is net sink */
	NEE = -1 * (NEP - fire);


	summary->NEP = NEP;
	summary->NPP = NPP;
	summary->NPPabove_nw = NPPabove_nw;
	summary->NPPbelow_nw = NPPbelow_nw;
	summary->NPPabove_w = NPPabove_w;
	summary->NPPbelow_w = NPPbelow_w;
	summary->NEE = NEE;
	summary->GPP = GPP;
	summary->MR = MR;
	summary->GR = GR;
	summary->HR = HR;
	summary->SR = SR;
	summary->TR = TR;
	summary->fire = fire;
	summary->cumNPP += NPP;
	summary->cumNPPabove_nw += NPPabove_nw;
	summary->cumNPPbelow_nw += NPPbelow_nw;
	summary->cumNPPabove_w += NPPabove_w;
	summary->cumNPPbelow_w += NPPbelow_w;
	summary->cumNEP += NEP;
	summary->cumNEE += NEE;
	summary->cumGPP += GPP;
	summary->cumMR += MR;
	summary->cumGR += GR;
	summary->cumHR += HR;
	summary->cumAR += AR;
	summary->cumTR += (MR + GR + HR);
	summary->cumSR += SR;
	summary->cumPET += wf->PET;
	summary->cumPETcanopy += wf->potETcanopy;
	summary->cumPEVPsurface += wf->potEVPandSUBLsurface;

	summary->cumN2Oflux    += summary->N2Oflux_total;
	summary->cumN2OfluxCeq += summary->N2OfluxCeq;
	summary->cumCH4flux    += summary->CH4flux_total;
	summary->cumCH4fluxCeq += summary->CH4fluxCeq;
	summary->cumEVPsurface += (wf->EVPsoilw + wf->EVPpondw + wf->SUBLsnoww);
	summary->cumETcanopy   += (wf->EVPcanopyw + wf->TRPsoilw_SUM);
	summary->cumET         += wf->ET;

	/************************************************************************************************************************************/
	/* 8. calculation litter fluxes and pools */

	summary->litdecomp = 0;
	for (layer = 0; layer < N_SOILLAYERS; layer++)
	{
		summary->litdecomp += cf->litr1c_to_soil1c[layer] + cf->litr2c_to_soil2c[layer] + cf->litr4c_to_soil3c[layer];
	}

	summary->litfire = cf->m_litr1c_to_fire_total + cf->m_litr2c_to_fire_total + cf->m_litr3c_to_fire_total + cf->m_litr4c_to_fire_total;

	/* aboveground littefall */
	summary->litfallc_above =
		cf->leafc_to_litr1c + cf->leafc_to_litr2c + cf->leafc_to_litr3c + cf->leafc_to_litr4c +
		cf->softstemc_to_litr1c + cf->softstemc_to_litr2c + cf->softstemc_to_litr3c + cf->softstemc_to_litr4c +
		cf->yieldc_to_litr1c + cf->yieldc_to_litr2c + cf->yieldc_to_litr3c + cf->yieldc_to_litr4c;

	/* belowground littefall*/
	summary->litfallc_below =
		cf->m_frootc_to_litr1c + cf->m_frootc_to_litr2c + cf->m_frootc_to_litr3c + cf->m_frootc_to_litr4c +
		cf->m_leafc_storage_to_litr1c + cf->m_frootc_storage_to_litr1c +
		cf->m_leafc_transfer_to_litr1c + cf->m_frootc_transfer_to_litr1c +
		cf->m_livestemc_storage_to_litr1c + cf->m_deadstemc_storage_to_litr1c +
		cf->m_livestemc_transfer_to_litr1c + cf->m_deadstemc_transfer_to_litr1c +
		cf->m_livecrootc_storage_to_litr1c + cf->m_deadcrootc_storage_to_litr1c +
		cf->m_livecrootc_transfer_to_litr1c + cf->m_deadcrootc_transfer_to_litr1c +
		cf->m_gresp_storage_to_litr1c + cf->m_gresp_transfer_to_litr1c +
		cf->frootc_to_litr1c + cf->frootc_to_litr2c + cf->frootc_to_litr3c + cf->frootc_to_litr4c +
		cf->m_softstemc_storage_to_litr1c + cf->m_softstemc_transfer_to_litr1c +
		cf->m_yieldc_storage_to_litr1c + cf->m_yieldc_transfer_to_litr1c;


	summary->litfallc = summary->litfallc_above + summary->litfallc_below;


	/************************************************************************************************************************************/
	/* 9. calculation of disturbance and senescence effect  */

	/* cut-down biomass and standing dead biome */
	STDB_to_litr = cf->STDBc_to_litr;
	CTDB_to_litr = cf->CTDBc_to_litr;

	summary->cumCplus_STDB += STDB_to_litr;
	summary->cumCplus_CTDB += CTDB_to_litr;

	/* management */
	Closs_THN_w = cf->livestemc_storage_to_THN + cf->livestemc_transfer_to_THN + cf->livestemc_to_THN +
		cf->deadstemc_storage_to_THN + cf->deadstemc_transfer_to_THN + cf->deadstemc_to_THN +
		cf->livecrootc_storage_to_THN + cf->livecrootc_transfer_to_THN + cf->livecrootc_to_THN +
		cf->deadcrootc_storage_to_THN + cf->deadcrootc_transfer_to_THN + cf->deadcrootc_to_THN;

	Closs_THN_nw = cf->leafc_storage_to_THN + cf->leafc_transfer_to_THN + cf->leafc_to_THN +
		cf->frootc_storage_to_THN + cf->frootc_transfer_to_THN + cf->frootc_to_THN +
		cf->yieldc_storage_to_THN + cf->yieldc_transfer_to_THN + cf->yieldc_to_THN +
		cf->gresp_transfer_to_THN + cf->gresp_storage_to_THN;


	Closs_MOW = cf->leafc_storage_to_MOW + cf->leafc_transfer_to_MOW + cf->leafc_to_MOW +
		cf->yieldc_storage_to_MOW + cf->yieldc_transfer_to_MOW + cf->yieldc_to_MOW +
		cf->softstemc_storage_to_MOW + cf->softstemc_transfer_to_MOW + cf->softstemc_to_MOW +
		cf->gresp_transfer_to_MOW + cf->gresp_storage_to_MOW;


	Closs_HRV = cf->leafc_storage_to_HRV + cf->leafc_transfer_to_HRV + cf->leafc_to_HRV +
		cf->yieldc_storage_to_HRV + cf->yieldc_transfer_to_HRV + cf->yieldc_to_HRV +
		cf->softstemc_storage_to_HRV + cf->softstemc_transfer_to_HRV + cf->softstemc_to_HRV +
		cf->gresp_transfer_to_HRV + cf->gresp_storage_to_HRV;


	yieldC_HRV = cf->yieldc_to_HRV;


	Closs_PLG = cf->leafc_storage_to_PLG - cf->leafc_transfer_to_PLG + cf->leafc_to_PLG +
		cf->yieldc_storage_to_PLG + cf->yieldc_transfer_to_PLG + cf->yieldc_to_PLG +
		cf->softstemc_storage_to_PLG - cf->softstemc_transfer_to_PLG - cf->softstemc_to_PLG +
		cf->frootc_storage_to_PLG + cf->frootc_transfer_to_PLG + cf->frootc_to_PLG +
		cf->gresp_transfer_to_PLG + cf->gresp_storage_to_PLG;


	Closs_GRZ = cf->leafc_storage_to_GRZ + cf->leafc_transfer_to_GRZ + cf->leafc_to_GRZ +
		cf->yieldc_storage_to_GRZ + cf->yieldc_transfer_to_GRZ + cf->yieldc_to_GRZ +
		cf->softstemc_storage_to_GRZ + cf->softstemc_transfer_to_GRZ + cf->softstemc_to_GRZ +
		cf->gresp_transfer_to_GRZ + cf->gresp_storage_to_GRZ;

	Cplus_GRZ = cf->GRZ_to_litr1c + cf->GRZ_to_litr2c + cf->GRZ_to_litr3c + cf->GRZ_to_litr4c;

	Cplus_FRZ = cf->FRZ_to_litr1c + cf->FRZ_to_litr2c + cf->FRZ_to_litr3c + cf->FRZ_to_litr4c;

	Cplus_PLT = cf->leafc_transfer_from_PLT + cf->frootc_transfer_from_PLT + cf->yieldc_transfer_from_PLT + cf->softstemc_transfer_from_PLT;

	Closs_PLT = cf->STDBc_leaf_to_PLT + cf->STDBc_froot_to_PLT + cf->STDBc_yield_to_PLT + cf->STDBc_softstem_to_PLT;


	Nplus_PLT = nf->leafn_transfer_from_PLT + nf->frootn_transfer_from_PLT + nf->yieldn_transfer_from_PLT + nf->softstemn_transfer_from_PLT;


	Nplus_GRZ = (nf->GRZ_to_litr1n + nf->GRZ_to_litr2n + nf->GRZ_to_litr3n + nf->GRZ_to_litr4n);
	Nplus_FRZ = (nf->FRZ_to_NH4 + nf->FRZ_to_NO3) +
		nf->FRZ_to_litr1n + nf->FRZ_to_litr2n + nf->FRZ_to_litr3n + nf->FRZ_to_litr4n;

	Nloss_PLT = nf->STDBn_leaf_to_PLT + nf->STDBn_froot_to_PLT + nf->STDBn_yield_to_PLT + nf->STDBn_softstem_to_PLT;

	summary->cumCloss_THN_w += Closs_THN_w;
	summary->cumCloss_THN_nw += Closs_THN_nw;
	summary->cumCloss_MOW += Closs_MOW;
	summary->cumCloss_HRV += Closs_HRV;
	summary->cumYieldC_HRV += yieldC_HRV;
	summary->cumCloss_PLG += Closs_PLG;
	summary->cumCloss_GRZ += Closs_GRZ;
	summary->cumCplus_GRZ += Cplus_GRZ;
	summary->cumCplus_PLT += Cplus_PLT;
	summary->cumCloss_PLT += Closs_PLT;
	summary->cumCplus_FRZ += Cplus_FRZ;
	summary->cumNplus_GRZ += Nplus_GRZ;
	summary->cumNplus_FRZ += Nplus_FRZ;
	summary->cumNplus_PLT += Nplus_PLT;
	summary->cumNloss_PLT += Nloss_PLT;
	summary->cumNplus_FRZ_org += nf->FRZ_to_litr1n + nf->FRZ_to_litr2n + nf->FRZ_to_litr3n + nf->FRZ_to_litr4n;
	summary->cumNplus_FRZ_NH4 += nf->FRZ_to_NH4 + nf->N2OfluxFRZ_NH4; //nf->FRZ_to_NH4: amount applied to the soil - volatilised part
	summary->cumNplus_FRZ_NO3 += nf->FRZ_to_NO3 + nf->N2OfluxFRZ_NO3;


	if (cs->vegCabove_HRV)
	{
		summary->harvestIndex = cs->yieldC_HRV / cs->vegCabove_HRV;
		summary->rootIndex = cs->frootC_HRV / (cs->frootC_HRV + cs->vegCabove_HRV);

	}
	else
	{
		summary->harvestIndex = 0;
		summary->rootIndex = 0;
	}

	if (epc->woody)
	{
		if (summary->frootc_LandD + summary->leafc_LandD + summary->yield_LandD + cs->livecrootc + cs->deadcrootc + cs->livestemc + cs->deadstemc)
			summary->belowground_ratio = (summary->frootc_LandD + cs->livecrootc + cs->deadcrootc) / (summary->frootc_LandD + summary->leafc_LandD + summary->yield_LandD + cs->livecrootc + cs->deadcrootc + cs->livestemc + cs->deadstemc);
		else
			summary->belowground_ratio = 0;
	}
	else
	{
		if (summary->frootc_LandD + summary->leafc_LandD + summary->yield_LandD + summary->softstemc_LandD)
			summary->belowground_ratio = summary->frootc_LandD / (summary->frootc_LandD + summary->leafc_LandD + summary->yield_LandD + summary->softstemc_LandD);
		else
			summary->belowground_ratio = 0;
	}

	/* senescence effect  */
	Closs_SNSC = cf->m_leafc_storage_to_SNSC + cf->m_leafc_transfer_to_SNSC + cf->m_leafc_to_SNSC +
		cf->m_yieldc_storage_to_SNSC + cf->m_yieldc_transfer_to_SNSC + cf->m_yieldc_to_SNSC +
		cf->m_softstemc_storage_to_SNSC + cf->m_softstemc_transfer_to_SNSC + cf->m_softstemc_to_SNSC +
		cf->m_frootc_storage_to_SNSC + cf->m_frootc_transfer_to_SNSC + cf->m_frootc_to_SNSC +
		cf->m_gresp_transfer_to_SNSC + cf->m_gresp_storage_to_SNSC +
		cf->HRV_frootc_to_SNSC + cf->HRV_softstemc_to_SNSC + cf->HRV_frootc_storage_to_SNSC + cf->HRV_frootc_transfer_to_SNSC +
		cf->HRV_softstemc_storage_to_SNSC + cf->HRV_softstemc_transfer_to_SNSC + cf->HRV_gresp_storage_to_SNSC + cf->HRV_gresp_transfer_to_SNSC;

	summary->cumCloss_SNSC += Closs_SNSC;

	/* NBP calculation: positive - mean net carbon gain to the system and negative - mean net carbon loss */
	disturb_loss = Closs_THN_w + Closs_THN_nw + Closs_MOW + Closs_HRV + Closs_PLG + Closs_GRZ;
	disturb_gain = Cplus_FRZ + Cplus_GRZ + Cplus_PLT;

	NBP = NEP + disturb_gain - disturb_loss;
	summary->NBP = NBP;
	summary->cumNBP += summary->NBP;

	/* lateral flux calculation */
	summary->cumCloss_MGM += disturb_loss;
	summary->cumCplus_MGM += disturb_gain;
	summary->Cflux_lateral = disturb_gain - disturb_loss;
	summary->cumCflux_lateral += summary->Cflux_lateral;

	/* NGB calculation: net greenhouse gas balance - NBP - N2O(Ceq) -CH(Ceq) */
	summary->NGB = summary->NBP - summary->N2OfluxCeq - summary->CH4fluxCeq;
	summary->cumNGB += summary->NGB;

	/* calculation of CN ratios */
	if (ns->litr1n_total > 0)
		epv->litr1_CNratio = cs->litr1c_total / ns->litr1n_total;
	else
		epv->litr1_CNratio = 0;

	if (ns->litr2n_total > 0)
		epv->litr2_CNratio = cs->litr2c_total / ns->litr2n_total;
	else
		epv->litr2_CNratio = 0;

	if (ns->litr3n_total > 0)
		epv->litr3_CNratio = cs->litr3c_total / ns->litr3n_total;
	else
		epv->litr3_CNratio = 0;

	if (ns->litr4n_total > 0)
		epv->litr4_CNratio = cs->litr4c_total / ns->litr4n_total;
	else
		epv->litr4_CNratio = 0;

	

	/************************************************************************************************************************************/
	/* Annual data for groundwater transported materail for total soil*/
	
	summary->cumGWmovchange_total += wf->GWmovchange_total;
	summary->cumGWtransp_total += wf->GWtransp_total;
	summary->cumGWdischarge += wf->GWdischarge;
	summary->cumGWrecharge += wf->GWrecharge;
	summary->cumGWevap += wf->GWevap;

	for (dm = 0; dm < N_DISSOLVMATER; dm++)
	{
		if (dm == 0)
			summary->cumNH4_GWchange_total += soilInfo->dismatGWmovchange[dm][layer] + soilInfo->dismatGWdecomp[dm][layer] + soilInfo->dismatGWecofunc[dm][layer] + soilInfo->dismatGWrecharge[dm] + soilInfo->dismatGWdischarge[dm];

		if (dm == 1)
			summary->cumNO3_GWchange_total += soilInfo->dismatGWmovchange[dm][layer] + soilInfo->dismatGWdecomp[dm][layer] + soilInfo->dismatGWecofunc[dm][layer] + soilInfo->dismatGWrecharge[dm] + soilInfo->dismatGWdischarge[dm];

		if (dm >= N_DISSOLVinorgN && dm < N_DISSOLVN)
			summary->cumSOILN_GWchange_total += soilInfo->dismatGWmovchange[dm][layer] + soilInfo->dismatGWdecomp[dm][layer] + soilInfo->dismatGWecofunc[dm][layer] + soilInfo->dismatGWrecharge[dm] + soilInfo->dismatGWdischarge[dm];

	}
	/************************************************************************************************************************************/
	/* Calculation balance data for unsaturated layers */

	/* initialization of storage variables */
	summary->NH4_unsat = 0;
	summary->NO3_unsat = 0;
	summary->orgN_unsat = vegN_above;
	summary->soilC_unsat = 0;
	summary->soilN_unsat = 0;


	/* management */
	summary->cumN2OfluxFRZ_NH4 += nf->N2OfluxFRZ_NH4;
	summary->cumN2OfluxFRZ_NO3 += nf->N2OfluxFRZ_NO3;
	summary->cumHRV_to_transpN += nf->HRV_to_transpN;
	summary->cumWflux_fromFRZ += wf->FRZ_to_soilw;
	summary->cumWflux_fromIRG += wf->IRG_to_prcp + wf->IRG_to_soilSurface + wf->IRG_to_soilw;


   /* noGW water src and snk */
	summary->cumWflux_fromPRCP += wf->prcp_to_canopyw + (wf->prcp_to_soilSurface - wf->IRG_to_prcp) + wf->prcp_to_snoww;
	summary->cumWflux_toET += wf->EVPcanopyw + wf->SUBLsnoww + wf->EVPsoilw + wf->TRPsoilw_SUM + wf->EVPpondw;
	summary->cumWflux_toRUNOFF += wf->pondw_to_runoff;


	/* GW water src and snk */
	summary->cumGWplus_unsat  += wf->GWevap + wf->GWtransp_zoneNORM + wf->GWtransp_zoneCAPIL + wf->GWdischarge;
	summary->cumGWminus_unsat += wf->GWrecharge;

	if (wf->GWmovchange_zoneNORM > 0)
		summary->cumGWplus_unsat += wf->GWmovchange_zoneNORM;
	else
		summary->cumGWminus_unsat += -1 * wf->GWmovchange_zoneNORM;

	if (wf->GWmovchange_zoneCAPIL > 0)
		summary->cumGWplus_unsat += wf->GWmovchange_zoneCAPIL;
	else
		summary->cumGWminus_unsat += -1 * wf->GWmovchange_zoneCAPIL;


	/* GW CandN src and snk: if no GW: unsat values are equal to total values */
	if (GWlayer != DATA_GAP)
	{
		/* above GWlayer */
		for (layer = 0; layer < GWlayer; layer++)
		{
			summary->NH4_unsat += ns->NH4[layer];
			summary->NO3_unsat += ns->NO3[layer];
			summary->orgN_unsat += vegN_below * epv->rootlengthProp[layer] + ns->litrN[layer] + ns->soilN[layer];
			summary->soilC_unsat += cs->soilC[layer];
			summary->soilN_unsat += ns->soilN[layer];

			summary->cumENV_to_NH4_unsat += nf->ndep_to_NH4[layer] + nf->nfix_to_NH4[layer];
			summary->cumENV_to_NO3_unsat += nf->ndep_to_NO3[layer];
			summary->cumNH4_to_soilSUM_unsat += nf->NH4_to_soilSUM[layer];
			summary->cumNO3_to_soilSUM_unsat += nf->NO3_to_soilSUM[layer];
			summary->cumNH4_to_nitrif_unsat += nf->NH4_to_nitrif[layer];
			summary->cumNO3_to_denitr_unsat += nf->NO3_to_denitr[layer];
			summary->cumNH4_to_npool_unsat += nf->NH4_to_npool[layer];
			summary->cumNO3_to_npool_unsat += nf->NO3_to_npool[layer];
			summary->cumN2OfluxNITRIF_unsat += nf->N2OfluxNITRIF[layer];
			

			for (dm = 0; dm < N_DISSOLVMATER; dm++)
			{
				/* NH4 */
				if (dm == 0)
				{ 
					if (soilInfo->dismatGWmovchange[dm][layer] > 0) 
						summary->cumNH4_GWplus_unsat += soilInfo->dismatGWmovchange[dm][layer];
					else
						summary->cumNH4_GWminus_unsat += -1 * soilInfo->dismatGWmovchange[dm][layer];

					if (soilInfo->dismatGWdecomp[dm][layer] > 0)
						summary->cumNH4_GWplus_unsat += soilInfo->dismatGWdecomp[dm][layer];
					else
						summary->cumNH4_GWminus_unsat += -1 * soilInfo->dismatGWdecomp[dm][layer];

					if (soilInfo->dismatGWecofunc[dm][layer] > 0)
						summary->cumNH4_GWplus_unsat += soilInfo->dismatGWecofunc[dm][layer];
					else
						summary->cumNH4_GWminus_unsat += -1 * soilInfo->dismatGWecofunc[dm][layer];
				}
				/* NO3 */
				if (dm == 1)
				{
		
					if (soilInfo->dismatGWmovchange[dm][layer] > 0)
						summary->cumNO3_GWplus_unsat += soilInfo->dismatGWmovchange[dm][layer];
					else
						summary->cumNO3_GWminus_unsat += -1 * soilInfo->dismatGWmovchange[dm][layer];

					if (soilInfo->dismatGWdecomp[dm][layer] > 0)
						summary->cumNO3_GWplus_unsat += soilInfo->dismatGWdecomp[dm][layer];
					else
						summary->cumNO3_GWminus_unsat += -1 * soilInfo->dismatGWdecomp[dm][layer];

					if (soilInfo->dismatGWecofunc[dm][layer] > 0)
						summary->cumNO3_GWplus_unsat += soilInfo->dismatGWecofunc[dm][layer];
					else
						summary->cumNO3_GWminus_unsat += -1 * soilInfo->dismatGWecofunc[dm][layer];
				}

				/* SOILN */
				if (dm >= N_DISSOLVinorgN && dm < N_DISSOLVN)
				{

					if (soilInfo->dismatGWmovchange[dm][layer] > 0)
						summary->cumSOILN_GWplus_unsat += soilInfo->dismatGWmovchange[dm][layer];
					else
						summary->cumSOILN_GWminus_unsat += -1 * soilInfo->dismatGWmovchange[dm][layer];

					if (soilInfo->dismatGWdecomp[dm][layer] > 0)
						summary->cumSOILN_GWplus_unsat += soilInfo->dismatGWdecomp[dm][layer];
					else
						summary->cumSOILN_GWminus_unsat += -1 * soilInfo->dismatGWdecomp[dm][layer];

					if (soilInfo->dismatGWecofunc[dm][layer] > 0)
						summary->cumSOILN_GWplus_unsat += soilInfo->dismatGWecofunc[dm][layer];
					else
						summary->cumSOILN_GWminus_unsat += -1 * soilInfo->dismatGWecofunc[dm][layer];
				}
	
			}
		}
		summary->cumNO3leach_unsat += soilInfo->dismatLeach[1][GWlayer - 1];
		summary->cumNO3leach_unsat = 0;
		/* in GWlayer */
		unSATprop = (sprop->dz_zoneNORM + sprop->dz_zoneCAPIL) / sitec->soillayer_thickness[GWlayer];

		summary->NH4_unsat += soilInfo->content_zoneNORM[0] + soilInfo->content_zoneCAPIL[0];
		summary->NO3_unsat += soilInfo->content_zoneNORM[1] + soilInfo->content_zoneCAPIL[1];
		summary->orgN_unsat += (vegN_below * epv->rootlengthProp[layer] + ns->litrN[layer]) * unSATprop;

		for (dm = N_DISSOLVinorgN; dm < N_DISSOLVMATER; dm++)
		{
			if (dm < N_DISSOLVN)
			{
				summary->orgN_unsat += soilInfo->content_zoneNORM[dm] + soilInfo->content_zoneCAPIL[dm];
				summary->soilN_unsat += soilInfo->content_zoneNORM[dm] + soilInfo->content_zoneCAPIL[dm];
			}
			else
				summary->soilC_unsat += soilInfo->content_zoneNORM[dm] + soilInfo->content_zoneCAPIL[dm];
		}
		
							summary->cumNH4_GWplus_unsat += soilInfo->dismatGWdischarge[dm];
					summary->cumNH4_GWminus_unsat += soilInfo->dismatGWrecharge[dm];


					summary->cumNO3_GWplus_unsat += soilInfo->dismatGWdischarge[dm];
					summary->cumNO3_GWminus_unsat += soilInfo->dismatGWrecharge[dm];

					summary->cumSOILN_GWplus_unsat += soilInfo->dismatGWdischarge[dm];
					summary->cumSOILN_GWminus_unsat += soilInfo->dismatGWrecharge[dm];

		summary->cumENV_to_NH4_unsat     += (nf->ndep_to_NH4[layer] + nf->nfix_to_NH4[layer]) * unSATprop;
		summary->cumENV_to_NO3_unsat     += (nf->ndep_to_NO3[layer])                          * unSATprop;
		summary->cumNH4_to_soilSUM_unsat += (nf->NH4_to_soilSUM[layer])                       * unSATprop;
		summary->cumNO3_to_soilSUM_unsat += (nf->NO3_to_soilSUM[layer])                       * unSATprop;
		summary->cumNH4_to_nitrif_unsat  += (nf->NH4_to_nitrif[layer])                        * unSATprop;
		summary->cumNO3_to_denitr_unsat  += (nf->NO3_to_denitr[layer])                        * unSATprop;
		summary->cumNH4_to_npool_unsat   += (nf->NH4_to_npool[layer])                         * unSATprop;
		summary->cumNO3_to_npool_unsat   += (nf->NO3_to_npool[layer])                         * unSATprop;
		summary->cumN2OfluxNITRIF_unsat  += (nf->N2OfluxNITRIF[layer])                        * unSATprop;

		for (dm = 0; dm < N_DISSOLVMATER; dm++)
		{
			// NH4 
			if (dm == 0)
			{
				summary->cumNH4_GWplus_unsat += soilInfo->dismatGWdischarge[dm];
				summary->cumNH4_GWminus_unsat += soilInfo->dismatGWrecharge[dm];

				if (soilInfo->dismatGWmovchange_NORM[dm] > 0)
					summary->cumNH4_GWplus_unsat += soilInfo->dismatGWmovchange_NORM[dm];
				else
					summary->cumNH4_GWminus_unsat += -1 * soilInfo->dismatGWmovchange_NORM[dm];

				if (soilInfo->dismatGWdecomp_NORM[dm] > 0)
					summary->cumNH4_GWplus_unsat += soilInfo->dismatGWdecomp_NORM[dm];
				else
					summary->cumNH4_GWminus_unsat += -1 * soilInfo->dismatGWdecomp_NORM[dm];

				if (soilInfo->dismatGWecofunc_NORM[dm] > 0)
					summary->cumNH4_GWplus_unsat += soilInfo->dismatGWecofunc_NORM[dm];
				else
					summary->cumNH4_GWminus_unsat += -1 * soilInfo->dismatGWecofunc_NORM[dm];

				if (soilInfo->dismatGWmovchange_CAPIL[dm] > 0)
					summary->cumNH4_GWplus_unsat += soilInfo->dismatGWmovchange_CAPIL[dm];
				else
					summary->cumNH4_GWminus_unsat += -1 * soilInfo->dismatGWmovchange_CAPIL[dm];

				if (soilInfo->dismatGWdecomp_CAPIL[dm] > 0)
					summary->cumNH4_GWplus_unsat += soilInfo->dismatGWdecomp_CAPIL[dm];
				else
					summary->cumNH4_GWminus_unsat += -1 * soilInfo->dismatGWdecomp_CAPIL[dm];

				if (soilInfo->dismatGWecofunc_CAPIL[dm] > 0)
					summary->cumNH4_GWplus_unsat += soilInfo->dismatGWecofunc_CAPIL[dm];
				else
					summary->cumNH4_GWminus_unsat += -1 * soilInfo->dismatGWecofunc_CAPIL[dm];
			}

			// NO3 
			if (dm == 1)
			{
				summary->cumNO3_GWplus_unsat += soilInfo->dismatGWdischarge[dm];
				summary->cumNO3_GWminus_unsat += soilInfo->dismatGWrecharge[dm];

				if (soilInfo->dismatGWmovchange_NORM[dm] > 0)
					summary->cumNO3_GWplus_unsat += soilInfo->dismatGWmovchange_NORM[dm];
				else
					summary->cumNO3_GWminus_unsat += -1 * soilInfo->dismatGWmovchange_NORM[dm];

				if (soilInfo->dismatGWdecomp_NORM[dm] > 0)
					summary->cumNO3_GWplus_unsat += soilInfo->dismatGWdecomp_NORM[dm];
				else
					summary->cumNO3_GWminus_unsat += -1 * soilInfo->dismatGWdecomp_NORM[dm];

				if (soilInfo->dismatGWecofunc_NORM[dm] > 0)
					summary->cumNO3_GWplus_unsat += soilInfo->dismatGWecofunc_NORM[dm];
				else
					summary->cumNO3_GWminus_unsat += -1 * soilInfo->dismatGWecofunc_NORM[dm];

				if (soilInfo->dismatGWmovchange_CAPIL[dm] > 0)
					summary->cumNO3_GWplus_unsat += soilInfo->dismatGWmovchange_CAPIL[dm];
				else
					summary->cumNO3_GWminus_unsat += -1 * soilInfo->dismatGWmovchange_CAPIL[dm];

				if (soilInfo->dismatGWdecomp_CAPIL[dm] > 0)
					summary->cumNO3_GWplus_unsat += soilInfo->dismatGWdecomp_CAPIL[dm];
				else
					summary->cumNO3_GWminus_unsat += -1 * soilInfo->dismatGWdecomp_CAPIL[dm];

				if (soilInfo->dismatGWecofunc_CAPIL[dm] > 0)
					summary->cumNO3_GWplus_unsat += soilInfo->dismatGWecofunc_CAPIL[dm];
				else
					summary->cumNO3_GWminus_unsat += -1 * soilInfo->dismatGWecofunc_CAPIL[dm];
			}
			// SOILN 
			if (dm >= N_DISSOLVinorgN && dm < N_DISSOLVN)
			{
				summary->cumSOILN_GWplus_unsat += soilInfo->dismatGWdischarge[dm];
				summary->cumSOILN_GWminus_unsat += soilInfo->dismatGWrecharge[dm];

				if (soilInfo->dismatGWmovchange_NORM[dm] > 0)
					summary->cumSOILN_GWplus_unsat += soilInfo->dismatGWmovchange_NORM[dm];
				else
					summary->cumSOILN_GWminus_unsat += -1 * soilInfo->dismatGWmovchange_NORM[dm];

				if (soilInfo->dismatGWdecomp_NORM[dm] > 0)
					summary->cumSOILN_GWplus_unsat += soilInfo->dismatGWdecomp_NORM[dm];
				else
					summary->cumSOILN_GWminus_unsat += -1 * soilInfo->dismatGWdecomp_NORM[dm];

				if (soilInfo->dismatGWecofunc_NORM[dm] > 0)
					summary->cumSOILN_GWplus_unsat += soilInfo->dismatGWecofunc_NORM[dm];
				else
					summary->cumSOILN_GWminus_unsat += -1 * soilInfo->dismatGWecofunc_NORM[dm];

				if (soilInfo->dismatGWmovchange_CAPIL[dm] > 0)
					summary->cumSOILN_GWplus_unsat += soilInfo->dismatGWmovchange_CAPIL[dm];
				else
					summary->cumSOILN_GWminus_unsat += -1 * soilInfo->dismatGWmovchange_CAPIL[dm];

				if (soilInfo->dismatGWdecomp_CAPIL[dm] > 0)
					summary->cumSOILN_GWplus_unsat += soilInfo->dismatGWdecomp_CAPIL[dm];
				else
					summary->cumSOILN_GWminus_unsat += -1 * soilInfo->dismatGWdecomp_CAPIL[dm];

				if (soilInfo->dismatGWecofunc_CAPIL[dm] > 0)
					summary->cumSOILN_GWplus_unsat += soilInfo->dismatGWecofunc_CAPIL[dm];
				else
					summary->cumSOILN_GWminus_unsat += -1 * soilInfo->dismatGWecofunc_CAPIL[dm];
			}

		}
	}
	else
	{
		summary->NH4_unsat = ns->NH4_total;
		summary->NO3_unsat = ns->NO3_total;
		summary->orgN_unsat = summary->orgN;
		summary->soilC_unsat = summary->soilC_total;
		summary->soilN_unsat = summary->soilN_total;

		for (layer = 0; layer < N_SOILLAYERS; layer++)
		{
			summary->cumENV_to_NH4_unsat += nf->ndep_to_NH4[layer] + nf->nfix_to_NH4[layer];
			summary->cumENV_to_NO3_unsat += nf->ndep_to_NO3[layer];
			summary->cumNH4_to_soilSUM_unsat += nf->NH4_to_soilSUM[layer];
			summary->cumNO3_to_soilSUM_unsat += nf->NO3_to_soilSUM[layer];
			summary->cumNH4_to_nitrif_unsat += nf->NH4_to_nitrif[layer];
			summary->cumNO3_to_denitr_unsat += nf->NO3_to_denitr[layer];
			summary->cumNH4_to_npool_unsat += nf->NH4_to_npool[layer];
			summary->cumNO3_to_npool_unsat += nf->NO3_to_npool[layer];
			summary->cumN2OfluxNITRIF_unsat += nf->N2OfluxNITRIF[layer];
		}
	}



	



	return(errorCode);
}
